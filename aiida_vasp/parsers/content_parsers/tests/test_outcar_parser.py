"""Test the OUTCAR parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
# pylint: disable=invalid-name

import numpy as np
import pytest

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


@pytest.mark.parametrize(['outcar_parser'], [(['disp_details', 'OUTCAR', {'quantities_to_parse': ['symmetries']}],)], indirect=True)
def test_parse_outcar_symmetry(outcar_parser, compare_symmetries):
    """Load a reference OUTCAR parser.

    We check that it parses and provides the correct content for the symmetry.

    """

    symmetry = outcar_parser.get_quantity('symmetries')
    assert set(symmetry) == set(compare_symmetries)


@pytest.mark.parametrize(['outcar_parser'], [(['disp_details', 'OUTCAR', {'quantities_to_parse': ['elastic_moduli']}],)], indirect=True)
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


@pytest.mark.parametrize(['outcar_parser'], [(['magnetization', 'OUTCAR', {
    'quantities_to_parse': ['magnetization', 'site_magnetization']
}],)],
                         indirect=True)
def test_parse_outcar_magnetizationr(fresh_aiida_env, outcar_parser):
    """Load a reference OUTCAR parser.

    We check that it parses and provides the correct content for the magnetization.

    """
    magnetization = outcar_parser.get_quantity('site_magnetization')
    test = {
        'sphere': {
            'x': {
                'site_moment': {
                    1: {
                        's': -0.014,
                        'p': -0.051,
                        'd': 1.687,
                        'tot': 1.621
                    },
                    2: {
                        's': -0.015,
                        'p': -0.052,
                        'd': 1.686,
                        'tot': 1.619
                    },
                    3: {
                        's': -0.014,
                        'p': -0.053,
                        'd': 1.708,
                        'tot': 1.64
                    },
                    4: {
                        's': -0.014,
                        'p': -0.053,
                        'd': 1.708,
                        'tot': 1.64
                    }
                },
                'total_magnetization': {
                    's': -0.057,
                    'p': -0.21,
                    'd': 6.788,
                    'tot': 6.521
                }
            },
            'y': {
                'site_moment': {},
                'total_magnetization': {}
            },
            'z': {
                'site_moment': {},
                'total_magnetization': {}
            }
        },
        'full_cell': np.asarray([6.4424922])
    }
    assert set(magnetization) == set(test)


@pytest.mark.parametrize(['outcar_parser'],
                         [(['magnetization_single', 'OUTCAR', {
                             'quantities_to_parse': ['magnetization', 'site_magnetization']
                         }],)],
                         indirect=True)
def test_parse_outcar_magnetization_single(fresh_aiida_env, outcar_parser):  # pylint: disable=invalid-name
    """Load a reference OUTCAR parser.

    We check that it parses and provides the correct content for the magnetization
    for a single atom in the cell.

    """
    magnetization = outcar_parser.get_quantity('site_magnetization')
    test = {
        'sphere': {
            'x': {
                'site_moment': {
                    1: {
                        's': -0.012,
                        'p': -0.043,
                        'd': 2.49,
                        'tot': 2.434
                    },
                },
                'total_magnetization': {
                    's': -0.012,
                    'p': -0.043,
                    'd': 2.49,
                    'tot': 2.434
                }
            },
            'y': {
                'site_moment': {},
                'total_magnetization': {}
            },
            'z': {
                'site_moment': {},
                'total_magnetization': {}
            }
        },
        'full_cell': np.asarray([2.4077611])
    }
    assert set(magnetization) == set(test)


@pytest.mark.parametrize('neb_outcar_parser', ['neb/01'], indirect=True)
def test_neb(fresh_aiida_env, neb_outcar_parser):
    """
    Test that the parameter node is a ParametersData instance.

    Should contain the symmetries and the elastic moduli.

    """
    data = neb_outcar_parser.get_quantity('neb_data')
    assert data['neb_converged']
    assert data['force_prep_real'] == pytest.approx(0.017467)
    assert data['energy_extrapolated'] == pytest.approx(-19.49550593)

    np.testing.assert_allclose(neb_outcar_parser.forces[0], np.array([0.008815, 0.005492, -0.000661]))
