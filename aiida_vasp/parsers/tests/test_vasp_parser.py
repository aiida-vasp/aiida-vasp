"""Unittests for VaspParser"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
# pylint: disable=protected-access,unused-variable,too-few-public-methods

import pytest

from aiida_vasp.io.parser import BaseFileParser
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_class, load_dbenv_if_not_loaded


class ExampleFileParser(BaseFileParser):
    """Example FileParser class for testing VaspParsers functionality."""

    PARSABLE_ITEMS = {
        'quantity_with_alternatives': {
            'inputs': [],
            'parsers': ['DUMMY'],
            'nodeName': 'structure',
            'prerequisites': [],
        },
        'quantity1': {
            'inputs': [],
            'parsers': ['CONTCAR'],
            'nodeName': 'structure',
            'is_alternative': 'quantity_with_alternatives',
            'prerequisites': []
        },
        'quantity2': {
            'inputs': [],
            'parsers': ['_scheduler-stdout.txt'],
            'nodeName': 'trajectory',
            'is_alternative': 'trajectory',
            'prerequisites': ['quantity1']
        },
        'quantity3': {
            'inputs': [],
            'parsers': ['CONTCAR'],
            'nodeName': '',
            'is_alternative': 'non_existing_quantity',
            'prerequisites': ['quantity_with_alternatives']
        },
    }

    def __init__(self, *args, **kwargs):
        super(ExampleFileParser, self).__init__(*args, **kwargs)
        self._parsable_items = ExampleFileParser.PARSABLE_ITEMS
        self._parsable_data = {}

    def _parse_file(self, inputs):
        from aiida.orm.data.parameter import ParameterData
        result = {}
        for quantity in ExampleFileParser.PARSABLE_ITEMS:
            result[quantity] = ParameterData(dict={})
        return result


class ExampleFileParser2(BaseFileParser):
    """Example class for testing non unique quantity identifiers."""

    PARSABLE_ITEMS = {
        'quantity1': {
            'inputs': [],
            'parsers': ['CONTCAR'],
            'nodeName': '',
            'is_alternative': 'quantity_with_alternatives',
            'prerequisites': []
        },
    }

    def __init__(self, *args, **kwargs):
        super(ExampleFileParser2, self).__init__(*args, **kwargs)
        self._parsable_items = ExampleFileParser2.PARSABLE_ITEMS
        self._parsable_data = {}

    def _parse_file(self, inputs):
        from aiida.orm.data.parameter import ParameterData
        result = {}
        for quantity in ExampleFileParser.PARSABLE_ITEMS:
            result[quantity] = ParameterData(dict={})
        return result


@pytest.fixture
def vasp_parser_with_test(vasp_nscf_and_ref, ref_retrieved_nscf):
    """Fixture providing a VaspParser instance coupled to a VaspCalculation."""
    from aiida.orm.data.parameter import ParameterData
    vasp_calc, _ = vasp_nscf_and_ref
    vasp_calc.use_settings(ParameterData(dict={'parser_settings': {'add_quantity_with_alternatives': True, 'add_quantity2': True}}))
    parser = vasp_calc.get_parserclass()(vasp_calc)
    parser.add_file_parser('_scheduler-stdout.txt', {'parser_class': ExampleFileParser, 'is_critical': False})
    success, outputs = parser.parse_with_retrieved({'retrieved': ref_retrieved_nscf})
    try:
        yield parser
    finally:
        parser = vasp_calc.get_parserclass()(vasp_calc)


@ONLY_ONE_CALC
def test_quantities_to_parse(vasp_nscf_and_ref, ref_retrieved_nscf):
    """Check if quantities are added to quantities to parse correctly."""
    from aiida.orm.data.parameter import ParameterData
    vasp_calc, _ = vasp_nscf_and_ref
    vasp_calc.use_settings(ParameterData(dict={'parser_settings': {'add_quantity_with_alternatives': True, 'add_quantity2': True}}))
    parser = vasp_calc.get_parserclass()(vasp_calc)

    parser.add_file_parser('_scheduler-stdout.txt', {'parser_class': ExampleFileParser, 'is_critical': False})

    parser.out_folder = parser.get_folder({'retrieved': ref_retrieved_nscf})
    parser._quantities.setup(parser._parsers, parser.out_folder)
    parser._parsers.setup(parser._settings)

    quantities_to_parse = parser._parsers.get_quantities_to_parse()
    assert 'quantity2' in quantities_to_parse
    assert 'quantity_with_alternatives' not in quantities_to_parse
    assert 'quantity1' in quantities_to_parse


@ONLY_ONE_CALC
def test_parsable_quantities(vasp_parser_with_test):
    """Check whether parsable quantities are set as intended."""
    parser = vasp_parser_with_test
    quantities = parser._quantities
    # Check whether all quantities from the added ExampleFileParser have been added.
    for quantity in ExampleFileParser.PARSABLE_ITEMS:
        assert quantities.get_by_name(quantity) is not None
    # Check whether quantities have been set up correctly.
    assert quantities.get_by_name('quantity1').has_files
    assert quantities.get_by_name('quantity1').is_parsable
    assert not quantities.get_by_name('quantity_with_alternatives').has_files
    assert quantities.get_by_name('quantity2').is_parsable
    assert not quantities.get_by_name('quantity3').is_parsable
    # check whether the additional non existing quantity has been added. This is for cases,
    # where a quantity is an alternative to another main quantity, which has not been loaded.
    assert quantities.get_by_name('non_existing_quantity') is not None


@ONLY_ONE_CALC
def test_node_link_names(vasp_parser_with_test):
    """Check whether an alternative quantity representing a node will be added with the correct linkname."""
    parser = vasp_parser_with_test
    assert 'quantity2' in parser._output_nodes
    # 'quantity2' is alternative to 'trajectory', which is not going to be parsed here.
    assert 'output_trajectory' in parser.new_nodes


@ONLY_ONE_CALC
def test_quantity_uniqeness(vasp_parser_with_test):
    """Make sure non-unique quantity identifiers are detected."""
    parser = vasp_parser_with_test
    # Add a second ExampleFileParser that defines a quantity with the same identifier as the first one.
    parser.add_file_parser('another_test_parser', {'parser_class': ExampleFileParser2, 'is_critical': False})
    with pytest.raises(RuntimeError) as excinfo:
        parser._quantities.setup(parser._parsers, parser.out_folder)
    assert 'quantity1' in str(excinfo.value)


def xml_path(folder):
    """Return the full path to the XML file."""
    return data_path(folder, 'vasprun.xml')


def xml_truncate(index, original, tmp):
    """Truncate vasprun.xml at the given line number and parse."""
    with open(original, 'r') as xmlfile:
        content = xmlfile.read().splitlines()
        truncated_content = '\n'.join(content[:-index or None])
    with open(tmp, 'w') as xmlfile:
        xmlfile.write(str(truncated_content))


@pytest.fixture(params=[0, 1])
def _parse_me(request, tmpdir):  # pylint disable=redefined-outer-name
    """
    Give the result of parsing a retrieved calculation (emulated).

    Returns a function which does:

    1. create a calculation with parser settings
    2. update the parser settings with the extra_settings
    3. create a parser with the calculation
    4. populate a fake retrieved folder and pass it to the parser
    5. return the result

    Note: This function is defined as protected to avoid pyling throwing
    redefined-outer-name (as it still does, even though it is disabled above).
    Pylint has some problems with pytest, there was a plugin, but maintanance is
    flaky. See discussions here:
    http://grokbase.com/t/python/pytest-dev/13bt9kz56y/fixtures-and-pylint-w0621


    """

    def parse(**extra_settings):
        """Run the parser using default settings updated with extra_settings."""
        load_dbenv_if_not_loaded()
        from aiida.orm import CalculationFactory, DataFactory
        from aiida_vasp.parsers.vasp import VaspParser
        calc = CalculationFactory('vasp.vasp')()
        settings_dict = {'parser_settings': {'add_bands': True, 'output_params': ['fermi_level']}}
        settings_dict.update(extra_settings)
        calc.use_settings(DataFactory('parameter')(dict=settings_dict))
        parser = VaspParser(calc=calc)
        retrieved = DataFactory('folder')()
        fldr = "basic"
        if "folder" in extra_settings:
            fldr = extra_settings["folder"]
        xml_file_path = xml_path(fldr)
        tmp_file_path = str(tmpdir.join('vasprun.xml'))
        #tmp_file_path = os.path.realpath(os.path.join(
        #    __file__, '../../../test_data/tmp/vasprun.xml'))
        xml_truncate(request.param, xml_file_path, tmp_file_path)
        retrieved.add_path(tmp_file_path, '')
        success, nodes = parser.parse_with_retrieved({'retrieved': retrieved})
        nodes = dict(nodes)
        return success, nodes

    return parse


def test_parser_nodes(_parse_me):
    """Test a few basic node items of the parser."""

    _, nodes = _parse_me(folder='basic')
    parameters = nodes['output_parameters']
    bands = nodes['output_bands']
    kpoints = nodes['output_kpoints']
    assert isinstance(parameters, get_data_class('parameter'))
    assert isinstance(bands, get_data_class('array.bands'))
    assert isinstance(kpoints, get_data_class('array.kpoints'))
    assert parameters.get_dict()['fermi_level'] == 5.96764939
