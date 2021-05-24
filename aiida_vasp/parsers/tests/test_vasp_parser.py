"""Unittests for VaspParser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
# pylint: disable=protected-access,unused-variable,too-few-public-methods, import-outside-toplevel

import os
import pytest
import numpy as np

from aiida.plugins import ParserFactory
from aiida.plugins import CalculationFactory
from aiida_vasp.parsers.file_parsers.parser import BaseFileParser
from aiida_vasp.parsers.vasp import NotificationComposer
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC, calc_with_retrieved
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_class


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
        from aiida.orm.nodes.data.dict import Dict
        result = {}
        for quantity in ExampleFileParser.PARSABLE_ITEMS:
            result[quantity] = Dict(dict={})
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
        from aiida.orm.nodes.data.dict import Dict
        result = {}
        for quantity in ExampleFileParser.PARSABLE_ITEMS:
            result[quantity] = Dict(dict={})
        return result


def _get_vasp_parser(calc_with_retrieved):
    """Return vasp parser before parsing"""
    settings_dict = {
        # 'ADDITIONAL_RETRIEVE_LIST': CalculationFactory('vasp.vasp')._ALWAYS_RETRIEVE_LIST,
        'parser_settings': {
            'add_custom': {
                'link_name': 'custom_node',
                'type': 'dict',
                'quantities': ['quantity2', 'quantity_with_alternatives']
            }
        }
    }
    file_path = str(os.path.abspath(os.path.dirname(__file__)) + '/../../test_data/basic_run')
    node = calc_with_retrieved(file_path, settings_dict)
    parser = ParserFactory('vasp.vasp')(node)
    return parser, file_path, node


@pytest.fixture
def vasp_parser_with_test(calc_with_retrieved):
    """Fixture providing a VaspParser instance coupled to a VaspCalculation."""
    parser, file_path, node = _get_vasp_parser(calc_with_retrieved)
    parser.add_parser_definition('_scheduler-stderr.txt', {'parser_class': ExampleFileParser, 'is_critical': False})
    success = parser.parse(retrieved_temporary_folder=file_path)
    try:
        yield parser
    finally:
        parser = ParserFactory('vasp.vasp')(node)


@pytest.fixture
def vasp_parser_without_parsing(calc_with_retrieved):
    parser, file_path, node = _get_vasp_parser(calc_with_retrieved)
    return parser, file_path


def test_add_parser_quantity_fail(vasp_parser_without_parsing):
    """add_parsable_quantity without file_name must fail"""
    parser, file_path = vasp_parser_without_parsing
    parser.add_parsable_quantity('quantity_with_alternatives', {
        'inputs': [],
        'prerequisites': [],
    })
    with pytest.raises(RuntimeError):
        parser.parse(retrieved_temporary_folder=file_path)


def test_add_parser_quantity(vasp_parser_without_parsing):
    """add_parsable_quantity with file_name succeeds."""
    parser, file_path = vasp_parser_without_parsing
    parser.add_parsable_quantity('quantity_with_alternatives', {'inputs': [], 'prerequisites': [], 'file_name': '_scheduler-stderr.txt'})
    parser.parse(retrieved_temporary_folder=file_path)
    quantities = parser._parsable_quantities
    assert 'quantity_with_alternatives' not in quantities.quantity_keys_to_parse
    assert 'quantity_with_alternatives' in quantities._missing_filenames
    assert quantities._missing_filenames['quantity_with_alternatives'] == '_scheduler-stderr.txt'


def test_add_parser_definition(vasp_parser_with_test):
    """Check if parser definition is passed to parser."""
    parser = vasp_parser_with_test
    parser_dict = parser._definitions.parser_definitions['_scheduler-stderr.txt']
    assert parser_dict['parser_class'] == ExampleFileParser


def test_parsable_quantities(vasp_parser_with_test):
    """Check whether parsable quantities are set as intended."""
    parser = vasp_parser_with_test
    quantity_keys_to_parse = parser._parsable_quantities.quantity_keys_to_parse
    missing_filenames = parser._parsable_quantities._missing_filenames
    quantity_keys = parser._parsable_quantities.quantity_keys_to_filenames.keys()
    # Check whether all quantities from the added ExampleFileParser have been added.
    for quantity_key in ExampleFileParser.PARSABLE_ITEMS:
        assert quantity_key in quantity_keys
    # Check whether quantities have been set up correctly.
    assert 'quantity1' not in missing_filenames
    assert 'quantity1' in quantity_keys_to_parse
    assert 'quantity2' not in quantity_keys_to_parse
    assert 'quantity3' not in quantity_keys_to_parse
    # check whether the additional non existing quantity has been added. This is for cases,
    # where a quantity is an alternative to another main quantity, which has not been loaded.
    assert 'non_existing_quantity' not in quantity_keys_to_parse


def test_quantity_uniqeness(vasp_parser_with_test):
    """Make sure non-unique quantity identifiers are detected."""
    parser = vasp_parser_with_test
    # Add a second ExampleFileParser that defines a quantity with the same identifier as the first one.
    parser.add_parser_definition('another_test_parser', {'parser_class': ExampleFileParser2, 'is_critical': False})
    with pytest.raises(RuntimeError) as excinfo:
        parser._parsable_quantities.setup(retrieved_filenames=parser._retrieved_content.keys(),
                                          parser_definitions=parser._definitions.parser_definitions,
                                          quantity_names_to_parse=parser._settings.quantity_names_to_parse)

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


def test_parser_nodes(request, calc_with_retrieved):
    """Test a few basic node items of the parser."""
    settings_dict = {'parser_settings': {'add_bands': True, 'add_kpoints': True, 'add_misc': ['fermi_level']}}

    file_path = str(request.fspath.join('..') + '../../../test_data/basic')

    node = calc_with_retrieved(file_path, settings_dict)

    parser_cls = ParserFactory('vasp.vasp')
    result, _ = parser_cls.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=file_path)

    misc = result['misc']
    bands = result['bands']
    kpoints = result['kpoints']

    assert isinstance(misc, get_data_class('dict'))
    assert isinstance(bands, get_data_class('array.bands'))
    assert isinstance(kpoints, get_data_class('array.kpoints'))
    assert misc.get_dict()['fermi_level'] == pytest.approx(5.96764939)


def test_parser_exception(request, calc_with_retrieved):
    """
    Test the handling of exceptions/missing quantities

    Here the eigenvalue section of the vasprun.xml and EIGNVAL files are deleted. However,
    we still expect the other propertie to be parsed correctly.
    """
    settings_dict = {
        'parser_settings': {
            'add_bands': True,
            'add_kpoints': True,
            'add_misc': ['total_energies', 'maximum_force', 'run_status', 'run_stats', 'notifications']
        }
    }

    file_path = str(request.fspath.join('..') + '../../../test_data/basic_run_ill_format')

    node = calc_with_retrieved(file_path, settings_dict)

    parser_cls = ParserFactory('vasp.vasp')
    result, output = parser_cls.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=file_path)

    assert output.is_finished
    assert output.exit_status == 1002

    misc = result['misc']
    assert isinstance(misc, get_data_class('dict'))
    assert misc.get_dict()['maximum_force'] == pytest.approx(0.0)
    assert misc.get_dict()['total_energies']['energy_extrapolated'] == pytest.approx(-36.09616894)

    assert misc['notifications'] == [{
        'name': "<class 'aiida_vasp.parsers.file_parsers.vasprun.VasprunParser'>",
        'message': 'the parser is not able to parse the occupancies quantity',
        'status': 1002,
    }]

    assert 'bands' not in result

    kpoints = result['kpoints']
    assert isinstance(kpoints, get_data_class('array.kpoints'))


def test_structure(request, calc_with_retrieved):
    """Test that the structure from vasprun and POSCAR is the same."""
    # turn of everything, except structure
    settings_dict = {
        'parser_settings': {
            'add_trajectory': False,
            'add_bands': False,
            'add_chgcar': False,
            'add_dos': False,
            'add_kpoints': False,
            'add_energies': False,
            'add_misc': False,
            'add_structure': True,
            'add_projectors': False,
            'add_born_charges': False,
            'add_dielectrics': False,
            'add_hessian': False,
            'add_dynmat': False,
            'add_wavecar': False,
            'add_site_magnetization': False,
        }
    }

    file_path = str(request.fspath.join('..') + '../../../test_data/basic')

    node = calc_with_retrieved(file_path, settings_dict)

    parser_cls = ParserFactory('vasp.vasp')
    result, _ = parser_cls.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=file_path)

    # First fetch structure from vasprun

    structure_vasprun = result['structure']
    assert isinstance(structure_vasprun, get_data_class('structure'))

    # Then from POSCAR/CONTCAR
    file_path = str(request.fspath.join('..') + '../../../test_data/basic_poscar')

    node = calc_with_retrieved(file_path, settings_dict)

    parser_cls = ParserFactory('vasp.vasp')
    result, _ = parser_cls.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=file_path)

    structure_poscar = result['structure']

    assert isinstance(structure_poscar, get_data_class('structure'))
    np.testing.assert_allclose(np.round(structure_vasprun.cell, 7), np.round(structure_poscar.cell, 7), rtol=0, atol=1e-8)
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
    np.testing.assert_allclose(positions_vasprun, positions_poscar, rtol=0, atol=1e-8)


def test_misc(request, calc_with_retrieved):
    """Test that it is possible to extract misc from both vasprun and OUTCAR."""
    # turn of everything, except misc
    settings_dict = {
        'parser_settings': {
            'add_trajectory': False,
            'add_bands': False,
            'add_chgcar': False,
            'add_dos': False,
            'add_kpoints': False,
            'add_energies': False,
            'add_misc': ['fermi_level', 'maximum_stress', 'maximum_force', 'total_energies', 'symmetries'],
            'add_structure': False,
            'add_projectors': False,
            'add_born_charges': False,
            'add_dielectrics': False,
            'add_hessian': False,
            'add_dynmat': False,
            'add_wavecar': False,
            'add_site_magnetization': False,
        }
    }

    file_path = str(request.fspath.join('..') + '../../../test_data/disp_details')

    node = calc_with_retrieved(file_path, settings_dict)

    parser_cls = ParserFactory('vasp.vasp')
    result, _ = parser_cls.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=file_path)

    misc = result['misc']
    assert isinstance(misc, get_data_class('dict'))
    data = misc.get_dict()
    # We already have a test to check if the quantities from the OUTCAR is correct, so
    # only perform rudimentary checks, and the content comming from the xml file.
    assert data['fermi_level'] == pytest.approx(6.17267267)
    assert data['maximum_stress'] == pytest.approx(42.96872956444064)
    assert data['maximum_force'] == pytest.approx(0.21326679)
    assert data['total_energies']['energy_extrapolated'] == pytest.approx(-10.823296)


@pytest.mark.parametrize(
    'config',
    [
        None,
        {
            'random_error': {
                'kind': 'ERROR',
                'regex': 'I AM A WELL DEFINED ERROR',
                'message': 'Okey, this error you do not want.',
                'suggestion': '',
                'location': 'STDOUT',
                'recover': False
            }
        },
        {
            'random_warning': {
                'kind': 'WARNING',
                'regex': 'I AM A WELL DEFINED WARNING',
                'message': 'Okey, this warning is nasty.',
                'suggestion': '',
                'location': 'STDOUT',
                'recover': False
            }
        }  # pylint: disable=too-many-statements, too-many-branches
    ])
@pytest.mark.parametrize('misc_input', [[], ['notifications']])
def test_stream(misc_input, config, request, calc_with_retrieved):
    """Test that the stream parser works and gets stored on a node."""
    file_path = str(request.fspath.join('..') + '../../../test_data/stdout/out')

    # turn of everything, except misc
    settings_dict = {
        'parser_settings': {
            'add_trajectory': False,
            'add_bands': False,
            'add_chgcar': False,
            'add_dos': False,
            'add_kpoints': False,
            'add_energies': False,
            'add_misc': misc_input,
            'add_structure': False,
            'add_projectors': False,
            'add_born_charges': False,
            'add_dielectrics': False,
            'add_hessian': False,
            'add_dynmat': False,
            'add_wavecar': False,
            'add_site_magnetization': False,
            'stream_config': config
        }
    }

    node = calc_with_retrieved(file_path, settings_dict)

    parser_cls = ParserFactory('vasp.vasp')
    result, _ = parser_cls.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=file_path)

    if misc_input == []:
        # Test empty misc specification, yields no misc output node
        with pytest.raises(KeyError) as error:
            misc = result['misc']
    else:
        misc = result['misc']
        misc_dict = misc.get_dict()
        if config is not None:
            if 'random_error' in config:
                assert len(misc_dict['notifications']) == 2
                assert misc_dict['notifications'][0]['name'] == 'ibzkpt'
                assert misc_dict['notifications'][0]['kind'] == 'ERROR'
                assert misc_dict['notifications'][0]['regex'] == 'internal error in subroutine IBZKPT'
                assert misc_dict['notifications'][1]['name'] == 'random_error'
                assert misc_dict['notifications'][1]['kind'] == 'ERROR'
                assert misc_dict['notifications'][1]['regex'] == 'I AM A WELL DEFINED ERROR'
            if 'random_warning' in config:
                assert len(misc_dict['notifications']) == 2
                assert misc_dict['notifications'][0]['name'] == 'ibzkpt'
                assert misc_dict['notifications'][0]['kind'] == 'ERROR'
                assert misc_dict['notifications'][0]['regex'] == 'internal error in subroutine IBZKPT'
                assert misc_dict['notifications'][1]['name'] == 'random_warning'
                assert misc_dict['notifications'][1]['kind'] == 'WARNING'
                assert misc_dict['notifications'][1]['regex'] == 'I AM A WELL DEFINED WARNING'
        else:
            assert len(misc_dict['notifications']) == 1
            assert misc_dict['notifications'][0]['name'] == 'ibzkpt'
            assert misc_dict['notifications'][0]['kind'] == 'ERROR'
            assert misc_dict['notifications'][0]['regex'] == 'internal error in subroutine IBZKPT'


def test_stream_history(request, calc_with_retrieved):
    """Test that the stream parser keeps history."""
    file_path = str(request.fspath.join('..') + '../../../test_data/stdout/out')

    # turn of everything, except misc
    settings_dict = {
        'parser_settings': {
            'add_trajectory': False,
            'add_bands': False,
            'add_chgcar': False,
            'add_dos': False,
            'add_kpoints': False,
            'add_energies': False,
            'add_misc': ['notifications'],
            'add_structure': False,
            'add_projectors': False,
            'add_born_charges': False,
            'add_dielectrics': False,
            'add_hessian': False,
            'add_dynmat': False,
            'add_wavecar': False,
            'add_site_magnetization': False,
            'stream_config': {
                'random_error': {
                    'kind': 'ERROR',
                    'regex': 'I AM A WELL DEFINED ERROR',
                    'message': 'Okey, this error you do not want.',
                    'suggestion': '',
                    'location': 'STDOUT',
                    'recover': False
                }
            },
            'stream_history': True
        }
    }

    node = calc_with_retrieved(file_path, settings_dict)

    parser_cls = ParserFactory('vasp.vasp')
    result, _ = parser_cls.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=file_path)

    misc = result['misc']
    misc_dict = misc.get_dict()
    assert len(misc_dict['notifications']) == 3
    assert misc_dict['notifications'][0]['name'] == 'ibzkpt'
    assert misc_dict['notifications'][0]['kind'] == 'ERROR'
    assert misc_dict['notifications'][0]['regex'] == 'internal error in subroutine IBZKPT'
    assert misc_dict['notifications'][1]['name'] == 'random_error'
    assert misc_dict['notifications'][1]['kind'] == 'ERROR'
    assert misc_dict['notifications'][1]['regex'] == 'I AM A WELL DEFINED ERROR'
    assert misc_dict['notifications'][2]['name'] == 'random_error'
    assert misc_dict['notifications'][2]['kind'] == 'ERROR'
    assert misc_dict['notifications'][2]['regex'] == 'I AM A WELL DEFINED ERROR'
    for item in misc_dict['notifications']:
        assert item['kind'] != 'WARNING'


def test_notification_composer(vasp_parser_without_parsing):
    """Test the NotificationComposer class"""
    parser, file_path = vasp_parser_without_parsing
    notifications = [{'name': 'edwav', 'kind': 'ERROR', 'message': 'Error in EDWAV'}]
    composer = NotificationComposer(notifications, {}, {'parameters': get_data_class('dict')(dict={'nelect': 10})}, parser.exit_codes)
    exit_code = composer.compose()
    assert exit_code.status == 703

    # BRMIX error but has NELECT defined in the input
    notifications = [{'name': 'brmix', 'kind': 'ERROR', 'message': 'Error in BRMIX'}]
    composer = NotificationComposer(notifications, {}, {'parameters': get_data_class('dict')(dict={'nelect': 10})}, parser.exit_codes)
    exit_code = composer.compose()
    assert exit_code is None

    # BRMIX error but no NELECT tag
    composer = NotificationComposer(notifications, {}, {'parameters': get_data_class('dict')(dict={})}, parser.exit_codes)
    exit_code = composer.compose()
    assert exit_code.status == 703
