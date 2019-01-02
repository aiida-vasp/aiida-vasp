"""Unittests for VaspParser"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
# pylint: disable=protected-access,unused-variable,too-few-public-methods

import pytest
import numpy as np

from aiida_vasp.io.parser import BaseFileParser
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_class, load_dbenv_if_not_loaded


class ExampleFileParser(BaseFileParser):
    """Example FileParser class for testing VaspParsers functionality."""

    PARSABLE_ITEMS = {
        'quantity1': {
            'inputs': [],
            'name': 'quantity_with_alternatives',
            'prerequisites': []
        },
        'quantity2': {
            'inputs': [],
            'name': 'trajectory',
            'prerequisites': ['quantity1']
        },
        'quantity3': {
            'inputs': [],
            'name': 'non_existing_quantity',
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
            'name': 'quantity_with_alternatives',
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
    vasp_calc.use_settings(
        ParameterData(
            dict={
                'parser_settings': {
                    'add_custom': {
                        'link_name': 'custom_node',
                        'type': 'parameter',
                        'quantities': ['quantity2', 'quantity_with_alternatives']
                    }
                }
            }))
    parser = vasp_calc.get_parserclass()(vasp_calc)
    parser.add_file_parser('_scheduler-stdout.txt', {'parser_class': ExampleFileParser, 'is_critical': False})
    parser.add_parsable_quantity(
        'quantity_with_alternatives',
        {
            'inputs': [],
            'prerequisites': [],
        },
    )
    success, outputs = parser.parse_with_retrieved({'retrieved': ref_retrieved_nscf})
    try:
        yield parser
    finally:
        parser = vasp_calc.get_parserclass()(vasp_calc)


@ONLY_ONE_CALC
def test_quantities_to_parse(vasp_parser_with_test):
    """Check if quantities are added to quantities to parse correctly."""
    parser = vasp_parser_with_test

    parser.quantities.setup()
    parser.parsers.setup()

    quantities_to_parse = parser.parsers.get_quantities_to_parse()
    assert 'quantity2' in quantities_to_parse
    assert 'quantity_with_alternatives' not in quantities_to_parse
    assert 'quantity1' in quantities_to_parse


@ONLY_ONE_CALC
def test_parsable_quantities(vasp_parser_with_test):
    """Check whether parsable quantities are set as intended."""
    parser = vasp_parser_with_test
    quantities = parser.quantities
    # Check whether all quantities from the added ExampleFileParser have been added.
    for quantity in ExampleFileParser.PARSABLE_ITEMS:
        assert quantities.get_by_name(quantity) is not None
    # Check whether quantities have been set up correctly.
    assert not quantities.get_by_name('quantity1').missing_files
    assert quantities.get_by_name('quantity1').is_parsable
    assert quantities.get_by_name('quantity_with_alternatives').missing_files
    assert quantities.get_by_name('quantity2').is_parsable
    assert not quantities.get_by_name('quantity3').is_parsable
    # check whether the additional non existing quantity has been added. This is for cases,
    # where a quantity is an alternative to another main quantity, which has not been loaded.
    assert quantities.get_by_name('non_existing_quantity') is not None


@ONLY_ONE_CALC
def test_quantity_uniqeness(vasp_parser_with_test):
    """Make sure non-unique quantity identifiers are detected."""
    parser = vasp_parser_with_test
    # Add a second ExampleFileParser that defines a quantity with the same identifier as the first one.
    parser.add_file_parser('another_test_parser', {'parser_class': ExampleFileParser2, 'is_critical': False})
    with pytest.raises(RuntimeError) as excinfo:
        parser.quantities.setup()
    assert 'quantity1' in str(excinfo.value)


def xml_path(folder):
    """Return the full path to the XML file."""
    return data_path(folder, 'vasprun.xml')


def poscar_path(folder):
    """Return the full path to the CONTCAR file."""
    return data_path(folder, 'CONTCAR')


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
        settings_dict = {'parser_settings': {'add_bands': True, 'add_kpoints': True, 'add_parameters': ['fermi_level']}}
        settings_dict.update(extra_settings)
        calc.use_settings(DataFactory('parameter')(dict=settings_dict))
        parser = VaspParser(calc=calc)
        retrieved = DataFactory('folder')()
        fldr = "basic"
        if "folder" in extra_settings:
            fldr = extra_settings["folder"]
        xml_file_path = xml_path(fldr)
        tmp_file_path = str(tmpdir.join('vasprun.xml'))
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


def test_structure(request):
    """Test that the structure from vasprun and POSCAR is the same."""

    load_dbenv_if_not_loaded()
    from aiida.orm import CalculationFactory, DataFactory
    from aiida_vasp.parsers.vasp import VaspParser
    calc = CalculationFactory('vasp.vasp')()
    # turn of everything, except structure
    settings_dict = {
        'parser_settings': {
            'add_trajectory': False,
            'add_bands': False,
            'add_chgcar': False,
            'add_dos': False,
            'add_kpoints': False,
            'add_energies': False,
            'add_parameters': False,
            'add_structure': True,
            'add_projectors': False,
            'add_born_charges': False,
            'add_dielectrics': False,
            'add_hessian': False,
            'add_dynmat': False,
            'add_wavecar': False,
            'file_parser_set': 'default'
        }
    }
    calc.use_settings(DataFactory('parameter')(dict=settings_dict))
    # First fetch structure from vasprun
    parser = VaspParser(calc=calc)
    retrieved = DataFactory('folder')()
    test_file_path = str(request.fspath.join('..') + '../../../test_data/basic/vasprun.xml')
    retrieved.add_path(test_file_path, '')
    success, nodes = parser.parse_with_retrieved({'retrieved': retrieved})
    nodes = dict(nodes)
    structure_vasprun = nodes['output_structure']
    assert isinstance(structure_vasprun, get_data_class('structure'))
    # Then from POSCAR/CONTCAR
    parser = VaspParser(calc=calc)
    retrieved = DataFactory('folder')()
    test_file_path = str(request.fspath.join('..') + '../../../test_data/basic_poscar/CONTCAR')
    retrieved.add_path(test_file_path, '')
    success, nodes = parser.parse_with_retrieved({'retrieved': retrieved})
    nodes = dict(nodes)
    structure_poscar = nodes['output_structure']
    assert isinstance(structure_poscar, get_data_class('structure'))
    assert np.array_equal(np.round(structure_vasprun.cell, 7), np.round(structure_poscar.cell, 7))
    positions_vasprun = []
    positions_poscar = []
    for site in structure_vasprun.sites:
        pos = np.round(np.asarray(site.position), 7)
        positions_vasprun.append(pos)
    for site in structure_poscar.sites:
        pos = np.round(np.asarray(site.position), 7)
        positions_poscar.append(pos)
    positions_vasprun = np.asarray(positions_vasprun)
    positions_poscar = np.asarray(positions_poscar)
    assert np.array_equal(positions_vasprun, positions_poscar)


def test_parameters(request):
    """Test that it is possible to extract parameters from both vasprun and OUTCAR."""

    load_dbenv_if_not_loaded()
    from aiida.orm import CalculationFactory, DataFactory
    from aiida_vasp.parsers.vasp import VaspParser
    calc = CalculationFactory('vasp.vasp')()
    # turn of everything, except parameters
    settings_dict = {
        'parser_settings': {
            'add_trajectory': False,
            'add_bands': False,
            'add_chgcar': False,
            'add_dos': False,
            'add_kpoints': False,
            'add_energies': False,
            'add_parameters': ['fermi_level', 'maximum_stress', 'maximum_force', 'total_energies', 'symmetries'],
            'add_structure': False,
            'add_projectors': False,
            'add_born_charges': False,
            'add_dielectrics': False,
            'add_hessian': False,
            'add_dynmat': False,
            'add_wavecar': False,
            'file_parser_set': 'default',
        }
    }
    calc.use_settings(DataFactory('parameter')(dict=settings_dict))
    parser = VaspParser(calc=calc)
    retrieved = DataFactory('folder')()
    vasprun_file_path = str(request.fspath.join('..') + '../../../test_data/disp_details/vasprun.xml')
    outcar_file_path = str(request.fspath.join('..') + '../../../test_data/disp_details/OUTCAR')
    retrieved.add_path(vasprun_file_path, '')
    retrieved.add_path(outcar_file_path, '')
    success, nodes = parser.parse_with_retrieved({'retrieved': retrieved})
    nodes = dict(nodes)
    parameters = nodes['output_parameters']
    assert isinstance(parameters, get_data_class('parameter'))
    data = parameters.get_dict()
    # We already have a test to check if the quantities from the OUTCAR is correct, so
    # only perform rudimentary checks, and the content comming from the xml file.
    assert data['fermi_level'] == 6.17267267
    assert data['maximum_stress'] == 42.96872956444064
    assert data['maximum_force'] == 0.21326679
    assert data['total_energies']['energy_no_entropy'] == -10.823296
