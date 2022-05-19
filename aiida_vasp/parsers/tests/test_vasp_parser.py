"""Unittests for VaspParser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
# pylint: disable=protected-access,unused-variable,too-few-public-methods, import-outside-toplevel
# pylint: disable=abstract-method

import os
from pathlib import Path
import pytest
import numpy as np

from aiida.plugins import ParserFactory
from aiida.plugins import CalculationFactory
from aiida.orm import load_node

from aiida_vasp.utils.aiida_utils import aiida_version, cmp_version
from aiida_vasp.parsers.content_parsers.base import BaseFileParser
from aiida_vasp.parsers.vasp import NotificationComposer
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC, calc_with_retrieved
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.parsers.content_parsers.kpoints import KpointsParser
from aiida_vasp.parsers.node_composer import NodeComposer


class ExampleFileParser(BaseFileParser):
    """Example object parser class for testing VaspParsers functionality."""

    PARSABLE_QUANTITIES = {
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
        super().__init__(*args, **kwargs)
        self._parsable_items = ExampleFileParser.PARSABLE_QUANTITIES
        self._parsable_data = {}

    def _parse_object(self, inputs):  # pylint: disable=no-self-use
        from aiida.orm.nodes.data.dict import Dict
        result = {}
        for quantity in ExampleFileParser.PARSABLE_QUANTITIES:
            result[quantity] = Dict(dict={})
        return result


class ExampleFileParser2(BaseFileParser):
    """Example class for testing non unique quantity identifiers."""

    PARSABLE_QUANTITIES = {
        'quantity1': {
            'inputs': [],
            'name': 'quantity_with_alternatives',
            'prerequisites': []
        },
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parsable_items = ExampleFileParser2.PARSABLE_QUANTITIES
        self._parsable_data = {}

    def _parse_object(self, inputs):  # pylint: disable=no-self-use
        from aiida.orm.nodes.data.dict import Dict
        result = {}
        for quantity in ExampleFileParser.PARSABLE_QUANTITIES:
            result[quantity] = Dict(dict={})
        return result


def _get_vasp_parser(calc_with_retrieved, request, settings_dict=None, relative_file_path=None):
    """Return vasp parser before parsing."""

    if settings_dict is None:
        _settings_dict = {
            # 'ADDITIONAL_RETRIEVE_LIST': CalculationFactory('vasp.vasp')._ALWAYS_RETRIEVE_LIST,
            'parser_settings': {
                'add_custom': {
                    'link_name': 'custom_node',
                    'type': 'dict',
                    'quantities': ['quantity2', 'quantity_with_alternatives']
                }
            }
        }
    else:
        _settings_dict = settings_dict
    if relative_file_path is None:
        _relative_file_path = '../../test_data/basic_run'
    else:
        _relative_file_path = relative_file_path

    # Path(request.fspath) will be replaced by request.node.path from pytest v7.
    file_path = str(Path(request.fspath).parent / _relative_file_path)
    node = calc_with_retrieved(file_path, _settings_dict)

    parser = ParserFactory('vasp.vasp')(node)
    return parser, file_path, node


@pytest.fixture
def vasp_parser_with_test(calc_with_retrieved, request):
    """Fixture providing a VaspParser instance coupled to a VaspCalculation."""
    parser, file_path, node = _get_vasp_parser(calc_with_retrieved, request)
    parser.add_parser_definition('_scheduler-stderr.txt', {'parser_class': ExampleFileParser, 'is_critical': False, 'mode': 'r'})
    success = parser.parse(retrieved_temporary_folder=file_path)
    try:
        yield parser
    finally:
        parser = ParserFactory('vasp.vasp')(node)


@pytest.fixture
def vasp_parser_without_parsing(calc_with_retrieved, request):
    parser, file_path, _ = _get_vasp_parser(calc_with_retrieved, request)
    return parser, file_path


def test_add_parser_quantity_fail(vasp_parser_without_parsing):
    """Add_parsable_quantity without name must fail"""
    parser, path = vasp_parser_without_parsing
    parser.add_parsable_quantity('quantity_with_alternatives', {
        'inputs': [],
        'prerequisites': [],
    })
    with pytest.raises(RuntimeError):
        parser.parse(retrieved_temporary_folder=path)


def test_add_parser_quantity(vasp_parser_without_parsing, aiida_caplog):
    """Add_parsable_quantity with name succeeds."""
    parser, path = vasp_parser_without_parsing
    parser.add_parsable_quantity('quantity_with_alternatives', {'inputs': [], 'prerequisites': [], 'name': '_scheduler-stderr.txt'})
    parser.parse(retrieved_temporary_folder=path)
    quantities = parser._parsable_quantities
    assert 'quantity_with_alternatives' not in quantities.quantity_keys_to_parse
    assert 'quantity_with_alternatives' in quantities._missing_content
    assert quantities._missing_content['quantity_with_alternatives'] == '_scheduler-stderr.txt'


def test_add_parser_definition(vasp_parser_with_test):
    """Check if parser definition is passed to parser."""
    parser = vasp_parser_with_test
    parser_dict = parser._definitions.parser_definitions['_scheduler-stderr.txt']
    assert parser_dict['parser_class'] == ExampleFileParser


def test_parsable_quantities(vasp_parser_with_test):
    """Check whether parsable quantities are set as intended."""
    parser = vasp_parser_with_test
    quantity_keys_to_parse = parser._parsable_quantities.quantity_keys_to_parse
    missing_content = parser._parsable_quantities._missing_content
    quantity_keys = parser._parsable_quantities.quantity_keys_to_content.keys()
    # Check whether all quantities from the added ExampleFileParser have been added.
    for quantity_key in ExampleFileParser.PARSABLE_QUANTITIES:
        assert quantity_key in quantity_keys
    # Check whether quantities have been set up correctly.
    assert 'quantity1' not in missing_content
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
    parser.add_parser_definition('another_test_parser', {'parser_class': ExampleFileParser2, 'is_critical': False, 'mode': 'r'})
    with pytest.raises(RuntimeError) as excinfo:
        parser._parsable_quantities.setup(retrieved_content=parser._retrieved_content.keys(),
                                          parser_definitions=parser._definitions.parser_definitions,
                                          quantity_names_to_parse=parser._settings.quantity_names_to_parse)

    assert 'quantity1' in str(excinfo.value)


def xml_path(folder):
    """Return the full path to the XML object."""
    return data_path(folder, 'vasprun.xml')


def poscar_path(folder):
    """Return the full path to the CONTCAR object."""
    return data_path(folder, 'CONTCAR')


def xml_truncate(index, original, tmp):
    """Truncate vasprun.xml at the given line number and parse."""
    with open(original, 'r', encoding='utf8') as handler:
        content = handler.read().splitlines()
        truncated_content = '\n'.join(content[:-index or None])
    with open(tmp, 'w', encoding='utf8') as handler:
        handler.write(str(truncated_content))


def test_parser_nodes(request, calc_with_retrieved):
    """Test a few basic node items of the parser."""
    settings_dict = {'parser_settings': {'add_bands': True, 'add_kpoints': True, 'add_misc': ['fermi_level']}}
    parser, file_path, _ = _get_vasp_parser(calc_with_retrieved,
                                            request,
                                            settings_dict=settings_dict,
                                            relative_file_path='../../test_data/basic')

    # The test data does not contain OUTCAR - make sure that is allowed
    parser._definitions.parser_definitions['OUTCAR']['is_critical'] = False
    parser.parse(retrieved_temporary_folder=file_path)
    result = parser.outputs

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
    we still expect the other properties to be parsed correctly.
    """

    settings_dict = {
        'parser_settings': {
            'add_bands': True,
            'add_kpoints': True,
            'add_misc': ['total_energies', 'maximum_force', 'run_status', 'run_stats', 'notifications']
        }
    }

    # Path(request.fspath) will be replaced by request.node.path from pytest v7.
    file_path = str(Path(request.fspath).parent / '../../test_data/basic_run_ill_format')
    node = calc_with_retrieved(file_path, settings_dict)
    parser_cls = ParserFactory('vasp.vasp')
    result, output = parser_cls.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=file_path)

    assert output.is_finished
    assert output.exit_status == 1002

    misc = result['misc']
    assert isinstance(misc, get_data_class('dict'))
    assert misc.get_dict()['maximum_force'] == pytest.approx(0.0)
    assert misc.get_dict()['total_energies']['energy_extrapolated'] == pytest.approx(-36.09616894)
    assert 'bands' not in result
    kpoints = result['kpoints']
    assert isinstance(kpoints, get_data_class('array.kpoints'))


@pytest.mark.xfail(aiida_version() < cmp_version('1.0.0a1'), reason='Element X only present in Aiida >= 1.x')
def test_parse_poscar_silly_read(fresh_aiida_env):
    """
    Parse (read) a reference POSCAR with silly elemental names.

    Using the PoscarParser and compare the result to a reference
    structure.

    """

    # Setup parser
    path = data_path('poscar', 'POSCARSILLY')
    parser = None
    with open(path, 'r', encoding='utf8') as handler:
        parser = PoscarParser(handler=handler)
    # Compose the node
    structure = parser.get_quantity('poscar-structure')
    result = NodeComposer.compose_structure('structure', {'structure': structure})
    names = result.get_site_kindnames()
    assert names == ['Hamburger', 'Pizza']
    symbols = result.get_symbols_set()
    assert symbols == set(['X', 'X'])


@pytest.mark.parametrize(['vasp_structure'], [('str-InAs',)], indirect=True)
def test_parse_poscar_silly_write(fresh_aiida_env, vasp_structure, tmpdir):
    """
    Parse (read, write) a reference POSCAR with silly elemental names.

    Using the PoscarParser and compare the result to a reference structure.

    """

    # Initialize the content parser
    parser = PoscarParser(data=vasp_structure)

    # Write silly POSCAR content to file
    temp_path = str(tmpdir.join('POSCAR'))
    parser.write(temp_path)

    # Read it again
    with open(temp_path, 'r', encoding='utf8') as handler:
        parser = PoscarParser(handler=handler)
    structure = parser.get_quantity('poscar-structure')
    result = NodeComposer.compose_structure('structure', {'structure': structure})

    # Compare
    names = result.get_site_kindnames()
    assert names == ['Hamburger', 'Pizza']
    symbols = result.get_symbols_set()
    assert symbols == set(['X', 'X'])


@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_parse_poscar_undercase(fresh_aiida_env, vasp_structure, tmpdir):
    """
    Parse a reference POSCAR.

    With potential elemental names using the PoscarParser and compare
    the result to a reference structure.

    """

    parser = PoscarParser(data=vasp_structure)
    result = parser.get_quantity('poscar-structure')
    names = result.get_site_kindnames()
    assert names == ['In', 'As', 'As', 'In_d', 'In_d', 'As']
    symbols = result.get_symbols_set()
    assert symbols == set(['As', 'In'])
    temp_path = str(tmpdir.join('POSCAR'))
    parser.write(temp_path)
    parser = None
    with open(temp_path, 'r', encoding='utf8') as handler:
        parser = PoscarParser(handler=handler)
    structure = parser.get_quantity('poscar-structure')
    result_reparse = NodeComposer.compose_structure('structure', {'structure': structure})
    names = result_reparse.get_site_kindnames()
    assert names == ['In', 'As', 'As', 'In_d', 'In_d', 'As']
    symbols = result_reparse.get_symbols_set()
    assert symbols == set(['As', 'In'])


def test_parse_kpoints(vasp_kpoints):
    """
    Parse a reference KPOINTS.

    Using the KpointsParser and compare the result to a reference
    KpointsData node.

    """

    kpoints, _ = vasp_kpoints

    try:
        _ = kpoints.get_attribute('mesh')
        path = data_path('kpoints', 'KPOINTS_mesh')
        method = 'get_kpoints_mesh'
        param = 'mesh'
    except AttributeError:
        pass

    try:
        _ = kpoints.get_attribute('array|kpoints')
        path = data_path('kpoints', 'KPOINTS_list')
        method = 'get_kpoints'
        param = 'list'
    except AttributeError:
        pass

    parser = None
    with open(path, 'r', encoding='utf8') as handler:
        parser = KpointsParser(handler=handler)
    kpts = parser.get_quantity('kpoints-kpoints')
    result = NodeComposer.compose_array_kpoints('array.kpoints', {'kpoints': kpts})
    if param == 'list':
        assert getattr(result, method)().all() == getattr(kpoints, method)().all()
    if param == 'mesh':
        assert getattr(result, method)() == getattr(kpoints, method)()


def test_structure(request, calc_with_retrieved):
    """Test that the structure from vasprun and POSCAR is the same."""
    # turn of everything, except structure
    settings_dict = {
        'parser_settings': {
            'add_trajectory': False,
            'add_bands': False,
            'add_charge_density': False,
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

    parser, file_path, _ = _get_vasp_parser(calc_with_retrieved,
                                            request,
                                            settings_dict=settings_dict,
                                            relative_file_path='../../test_data/basic')

    # The test data does not contain OUTCAR - make sure that is allowed
    parser._definitions.parser_definitions['OUTCAR']['is_critical'] = False
    parser.parse(retrieved_temporary_folder=file_path)
    result = parser.outputs

    # First fetch structure from vasprun
    structure_vasprun = result['structure']
    assert isinstance(structure_vasprun, get_data_class('structure'))

    # Then from POSCAR/CONTCAR
    parser, file_path, _ = _get_vasp_parser(calc_with_retrieved,
                                            request,
                                            settings_dict=settings_dict,
                                            relative_file_path='../../test_data/basic_poscar')
    # The test data does not contain OUTCAR - make sure that is allowed
    parser._definitions.parser_definitions['OUTCAR']['is_critical'] = False
    parser._definitions.parser_definitions['vasprun.xml']['is_critical'] = False
    parser.parse(retrieved_temporary_folder=file_path)
    result = parser.outputs

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
            'add_charge_density': False,
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

    # Path(request.fspath) will be replaced by request.node.path from pytest v7.
    file_path = str(Path(request.fspath).parent / '../../test_data/disp_details')
    node = calc_with_retrieved(file_path, settings_dict)
    parser_cls = ParserFactory('vasp.vasp')
    result, _ = parser_cls.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=file_path)

    misc = result['misc']
    assert isinstance(misc, get_data_class('dict'))
    data = misc.get_dict()
    # We already have a test to check if the quantities from the OUTCAR is correct, so
    # only perform rudimentary checks, and the content comming from the xml object.
    assert data['fermi_level'] == pytest.approx(6.17267267)
    assert data['maximum_stress'] == pytest.approx(42.96872956444064)
    assert data['maximum_force'] == pytest.approx(0.21326679)
    assert data['total_energies']['energy_extrapolated'] == pytest.approx(-10.823296)


def test_custom_outputs(request, calc_with_retrieved):
    """Test custom_outputs by fermi_level."""
    parser_settings = {
        'add_custom_outputs': {
            'type': 'float',
            'quantities': ['fermi_level',],
            'link_name': 'custom_outputs.fermi_level',
        },
        'add_custom_outputs2': {
            'type': 'float',
            'quantities': ['fermi_level',],
            'link_name': 'fermi_level3',
        },
        'add_fermi_level2': {
            'type': 'float',
            'quantities': ['fermi_level',],
        },
        'add_fermi_level': {
            'type': 'float'
        }
    }
    for parser_setting in parser_settings:
        parser, file_path, _ = _get_vasp_parser(calc_with_retrieved, request, settings_dict={'parser_settings': parser_settings})
        parser.parse(retrieved_temporary_folder=file_path)
        assert 'custom_outputs.fermi_level' in parser.outputs
        assert 'fermi_level3' in parser.outputs
        assert 'fermi_level2' in parser.outputs
        assert 'fermi_level' in parser.outputs
        assert parser.outputs['custom_outputs.fermi_level'] == pytest.approx(4.29634683)


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
    # turn of everything, except misc
    settings_dict = {
        'parser_settings': {
            'add_trajectory': False,
            'add_bands': False,
            'add_charge_density': False,
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

    parser, file_path, _ = _get_vasp_parser(calc_with_retrieved,
                                            request,
                                            settings_dict=settings_dict,
                                            relative_file_path='../../test_data/stdout/out')
    parser._definitions.parser_definitions['OUTCAR']['is_critical'] = False
    parser._definitions.parser_definitions['vasprun.xml']['is_critical'] = False
    parser.parse(retrieved_temporary_folder=file_path)
    result = parser.outputs
    assert result is not None

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
    # turn of everything, except misc
    settings_dict = {
        'parser_settings': {
            'add_trajectory': False,
            'add_bands': False,
            'add_charge_density': False,
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

    parser, file_path, _ = _get_vasp_parser(calc_with_retrieved,
                                            request,
                                            settings_dict=settings_dict,
                                            relative_file_path='../../test_data/stdout/out')

    # The test data does not contain OUTCAR - make sure that is allowed
    parser._definitions.parser_definitions['OUTCAR']['is_critical'] = False
    parser._definitions.parser_definitions['vasprun.xml']['is_critical'] = False
    parser.parse(retrieved_temporary_folder=file_path)
    result = parser.outputs

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
    parser, path = vasp_parser_without_parsing
    notifications = [{'name': 'edwav', 'kind': 'ERROR', 'message': 'Error in EDWAV'}]
    composer = NotificationComposer(notifications, {}, {'parameters': get_data_class('dict')(dict={
        'nelect': 10
    })},
                                    parser.exit_codes,
                                    parser_settings=parser._settings)
    exit_code = composer.compose()
    assert exit_code.status == 703

    # BRMIX error but has NELECT defined in the input
    notifications = [{'name': 'brmix', 'kind': 'ERROR', 'message': 'Error in BRMIX'}]
    composer = NotificationComposer(notifications, {}, {'parameters': get_data_class('dict')(dict={
        'nelect': 10
    })},
                                    parser.exit_codes,
                                    parser_settings=parser._settings)
    exit_code = composer.compose()
    assert exit_code is None

    # BRMIX error but no NELECT tag
    composer = NotificationComposer(notifications, {}, {'parameters': get_data_class('dict')(dict={})},
                                    parser.exit_codes,
                                    parser_settings=parser._settings)
    exit_code = composer.compose()
    assert exit_code.status == 703


def test_critical_object_missing(calc_with_retrieved, request):
    """Test raising return code to indicate that one or more critical objects are missing"""
    # Here we add a fictional file as critical
    parser, file_path, node = _get_vasp_parser(calc_with_retrieved, request)
    parser.add_parser_definition('some-critical-file.txt', {'parser_class': ExampleFileParser, 'is_critical': True, 'mode': 'r'})
    parser.add_parsable_quantity('quantity_with_alternatives', {'inputs': [], 'prerequisites': [], 'file_name': '_scheduler-stderr.txt'})
    success = parser.parse(retrieved_temporary_folder=file_path)
    assert success == parser.exit_codes.ERROR_CRITICAL_MISSING_OBJECT

    # Here we remove a file that is marked as critical by default
    parser, file_path, node = _get_vasp_parser(calc_with_retrieved, request)
    # Delete the retrieved OUTCAR file and instantiate the parser. In order to
    # delete from a stored node we need to access the backend, which allows such deletion.
    # Also, make sure we use the same variable here, otherwise we would e.g. issue the
    # update_repository_metadata on a different repository entry in memory.
    retrieved = node.outputs.retrieved
    retrieved.base.repository._repository.delete_object('OUTCAR')
    retrieved.base.repository._update_repository_metadata()
    parser = ParserFactory('vasp.vasp')(node)
    # No temporary folder is passed - parse everything for the permanent storage
    success = parser.parse()
    assert success == parser.exit_codes.ERROR_CRITICAL_MISSING_OBJECT
