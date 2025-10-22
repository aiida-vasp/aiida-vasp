"""Unittests for VtstNebVaspParser."""

import pathlib

import numpy as np
import pytest
from aiida.plugins import ParserFactory

cwd = pathlib.Path(__file__).parent


def _get_neb_vasp_parser(neb_calc_with_retrieved):
    """Return vasp parser before parsing"""
    settings_dict = {
        # 'ADDITIONAL_RETRIEVE_LIST': CalculationFactory('vasp.vasp')._ALWAYS_RETRIEVE_LIST,
        'parser_settings': {}
    }
    file_path = cwd / '..' / 'test_data/neb'
    node = neb_calc_with_retrieved(file_path, settings_dict, 3)
    parser = ParserFactory('vasp.neb')(node)
    return parser, file_path, node


@pytest.fixture
def neb_parser_with_test(neb_calc_with_retrieved):
    """Fixture providing a VaspParser instance coupled to a VaspCalculation."""
    parser, _, node = _get_neb_vasp_parser(neb_calc_with_retrieved)
    try:
        yield parser
    finally:
        parser = ParserFactory('vasp.vasp')(node)


def test_neb_parser(neb_parser_with_test):
    """
    Test the neb parser
    """
    parser = neb_parser_with_test
    parser.parse()
    assert 'misc' in neb_parser_with_test.outputs
    misc_dict = parser.outputs.misc.get_dict()

    assert misc_dict['neb_data']['01']['neb_converged']
    assert misc_dict['total_energies']['01']['energy_free'] == -19.49164066
    assert '02' in misc_dict['total_energies']

    assert misc_dict['neb_data']['03']

    # Check that notifications exists
    assert 'notifications' in misc_dict

    # Make sure structures are parsed as well
    assert 'structure.image_01' in parser.outputs

    # Check that the forces is parsed
    forces = misc_dict['forces']['01']
    assert np.all(forces[0] == np.array([0.008815, 0.005492, -0.000661]))

    assert np.array(forces).shape == (4, 3)
