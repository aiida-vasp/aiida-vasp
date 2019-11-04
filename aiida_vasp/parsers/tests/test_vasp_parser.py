"""Unittests for VaspParser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
# pylint: disable=protected-access,unused-variable,too-few-public-methods

import os
import pytest
import numpy as np

from aiida_vasp.parsers.file_parsers.parser import BaseFileParser
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


@pytest.fixture
def vasp_parser_with_test(calc_with_retrieved):
    """Fixture providing a VaspParser instance coupled to a VaspCalculation."""
    from aiida.plugins import ParserFactory
    from aiida.plugins import CalculationFactory

    settings_dict = {
        #'ADDITIONAL_RETRIEVE_LIST': CalculationFactory('vasp.vasp')._ALWAYS_RETRIEVE_TEMPORARY_LIST,
        'parser_settings': {
            'add_custom': {
                'link_name': 'custom_node',
                'type': 'dict',
                'quantities': ['quantity2', 'quantity_with_alternatives']
            }
        }
    }

    file_path = str(os.path.abspath(os.path.dirname(__file__)) + '/../../test_data/test_relax_wc/out')

    node = calc_with_retrieved(file_path, settings_dict)

    parser = ParserFactory('vasp.vasp')(node)
    parser.add_file_parser('_scheduler-stderr.txt', {'parser_class': ExampleFileParser, 'is_critical': False})
    parser.add_parsable_quantity(
        'quantity_with_alternatives',
        {
            'inputs': [],
            'prerequisites': [],
        },
    )
    success = parser.parse(retrieved_temporary_folder=file_path)
    try:
        yield parser
    finally:
        parser = ParserFactory('vasp.vasp')(node)


def test_quantities_to_parse(vasp_parser_with_test):
    """Check if quantities are added to quantities to parse correctly."""
    parser = vasp_parser_with_test

    parser.quantities.setup()
    parser.parsers.setup()

    quantities_to_parse = parser.parsers.get_quantities_to_parse()
    assert 'quantity2' in quantities_to_parse
    assert 'quantity_with_alternatives' not in quantities_to_parse
    assert 'quantity1' in quantities_to_parse


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


def test_parser_nodes(request, calc_with_retrieved):
    """Test a few basic node items of the parser."""
    from aiida.plugins import ParserFactory

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
    assert misc.get_dict()['fermi_level'] == 5.96764939


def test_structure(request, calc_with_retrieved):
    """Test that the structure from vasprun and POSCAR is the same."""
    from aiida.plugins import ParserFactory

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
            'file_parser_set': 'default'
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


def test_misc(request, calc_with_retrieved):
    """Test that it is possible to extract misc from both vasprun and OUTCAR."""
    from aiida.plugins import ParserFactory

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
            'file_parser_set': 'default',
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
    assert data['fermi_level'] == 6.17267267
    assert data['maximum_stress'] == 42.96872956444064
    assert data['maximum_force'] == 0.21326679
    assert data['total_energies']['energy_no_entropy'] == -10.823296
