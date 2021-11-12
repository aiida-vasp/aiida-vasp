"""Test the OUTCAR parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy as np

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path


@pytest.mark.parametrize('outcar_parser', ['outcar'], indirect=True)
def test_parse_outcar_bogus_quantity(outcar_parser):
    """Load a reference OUTCAR parser.

    We check that the get quantity returns None for bogus keys.

    """

    assert outcar_parser.get_quantity('bogusquantity1234') is None


@pytest.mark.parametrize('outcar_parser', ['outcar'], indirect=True)
def test_parse_outcar_no_quantity(outcar_parser):
    """Load a reference OUTCAR parser.

    We check that the get quantity returns None if that quantity is not found
    in the presented OUTCAR.

    """
    # A standard OUTCAR should not contain parsable magnetization.
    assert outcar_parser.get_quantity('magnetization') is None


@pytest.mark.parametrize('outcar_parser', ['disp_details'], indirect=True)
def test_parse_outcar_symmetry(outcar_parser):
    """Load a reference OUTCAR parser.

    We check that it parses and provides the correct content for the symmetry.

    """

    symmetry = outcar_parser.get_quantity('symmetries')
    # test symmetries
    test = {
        'symmetrized_cell_type': {
            'static': [
                'face centered cubic supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'face centered cubic supercell.', 'face centered cubic supercell.',
                'face centered cubic supercell.'
            ],
            'dynamic': [
                'face centered cubic supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'face centered cubic supercell.', 'face centered cubic supercell.',
                'face centered cubic supercell.'
            ]
        },
        'original_cell_type': {
            'static': [
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell'
            ],
            'dynamic': [
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell'
            ]
        },
        'num_space_group_operations': {
            'static': [48, 16, 16, 16, 16, 16, 16, 4, 4, 4, 4, 4, 4, 8, 8, 48],
            'dynamic': [48, 16, 16, 16, 16, 16, 16, 4, 4, 4, 4, 4, 4, 8, 8, 48]
        }
    }
    assert set(symmetry) == set(test)


@pytest.mark.parametrize('outcar_parser', ['disp_details'], indirect=True)
def test_parse_outcar_elastic_moduli(outcar_parser):
    """Load a reference OUTCAR parser.

    We check that it parses and provides the correct content for the elastic moduli.

    """

    moduli = outcar_parser.get_quantity('elastic_moduli')
    test = np.array([1674.5786, 704.739, 704.739, -0.0, 0.0, 0.0])
    np.testing.assert_allclose(moduli['symmetrized'][0], test)
    test = np.array([0.0, 0.0, 0.0, -0.0, -0.0, 1122.6622])
    np.testing.assert_allclose(moduli['symmetrized'][5], test)
    test = np.array([705.0238, 1674.8491, 705.0238, -0.0, -0.0, 0.0])
    np.testing.assert_allclose(moduli['non_symmetrized'][1], test)
    test = np.array([-0.0078, -0.0495, 0.0147, 0.0, 1123.0829, -0.0])
    np.testing.assert_allclose(moduli['non_symmetrized'][4], test)
    test = np.array([704.739, 704.739, 1674.5786, -0.0, -0.0, 0.0])
    np.testing.assert_allclose(moduli['total'][2], test)
    test = np.array([-0.0, -0.0, -0.0, 775.8054, 0.0, -0.0])
    np.testing.assert_allclose(moduli['total'][3], test)


@pytest.mark.parametrize('outcar_parser', ['disp_details'], indirect=True)
def test_parse_outcar_status(outcar_parser):
    """Load a reference OUTCAR parser.

    We check that it parses and provides the correct content for a simple status.

    """

    status = outcar_parser.get_quantity('run_status')
    assert status['last_iteration_index'] == [15, 5]
    assert status['finished']
    assert status['ionic_converged']
    assert status['electronic_converged']
    assert status['nelm'] == 60
    assert status['nsw'] == 61


_TEST_DATA = [
    (['outcar_extras', 'OUTCAR.converged'], [True, True, True, False, False]),
    (['outcar_extras', 'OUTCAR.nelm-breach-consistent'], [True, False, False, True, True]),
    (['outcar_extras', 'OUTCAR.nelm-breach-partial'], [True, False, True, False, True]),
    (['outcar_extras', 'OUTCAR.unfinished'], [False, False, False, False, False]),
    (['outcar_extras', 'OUTCAR.not-converged'], [True, False, True, False, False]),
]


@pytest.mark.parametrize('outcar_parser,expected', _TEST_DATA, indirect=['outcar_parser'])
def test_parse_outcar_status_extended(outcar_parser, expected):
    """Load a reference OUTCAR parser.

    We check that it parses and provides the correct content for the status for multiple OUTCARs.

    """

    status = outcar_parser.get_quantity('run_status')
    assert status['finished'] is expected[0]
    assert status['ionic_converged'] is expected[1]
    assert status['electronic_converged'] is expected[2]
    assert status['consistent_nelm_breach'] is expected[3]
    assert status['contains_nelm_breach'] is expected[4]
