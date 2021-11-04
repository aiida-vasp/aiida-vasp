"""Test the DOSCAR parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy as np

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path


@pytest.mark.parametrize(['doscar_parser'], [('doscar',)], indirect=True)
def test_parse_doscar(fresh_aiida_env, doscar_parser):
    """Load a reference DOSCAR parser.

    We check that the it parses and provides the correct content.

    """
    result = doscar_parser.get_quantity('doscar-dos')
    compare_doscar_content(result)


@pytest.mark.parametrize(['doscar_parser'], [(['doscar', 'DOSCAR.nopdos'], )], indirect=True)
def test_parse_doscar_nodos(fresh_aiida_env, doscar_parser):
    """Load a reference DOSCAR parser with no partial density of states.

    We check that the it parses and provides the correct content.

    """
    result = doscar_parser.get_quantity('doscar-dos')
    compare_doscar_content(result, pdos=False)

@pytest.mark.parametrize(['doscar_parser'], [(['doscar', 'DOSCAR.spin'], )], indirect=True)
def test_parse_doscar_spin(fresh_aiida_env, doscar_parser):
    """Parse a reference DOSCAR with the DosParser and compare the result to a reference."""
    result = doscar_parser.get_quantity('doscar-dos')

    result_dos = result['tdos']
    assert len(result_dos.dtype) == 3
    assert result_dos['total'].shape == (301, 2)

    result_dos = result['pdos']
    assert len(result_dos.dtype) == 10
    assert result_dos[0]['px'].shape == (301, 2)

@pytest.mark.parametrize(['doscar_parser'], [(['doscar', 'DOSCAR.ncl'], )], indirect=True)
def test_parse_doscar_ncl(fresh_aiida_env, doscar_parser):
    """parse a reference DOSCAR with the dosparser and compare the result to a reference."""
    result = doscar_parser.get_quantity('doscar-dos')

    result_dos = result['tdos']
    assert len(result_dos.dtype) == 3
    assert result_dos['total'].shape == (301,)

    result_dos = result['pdos']
    assert len(result_dos.dtype) == 10
    assert result_dos[0]['px'].shape == (301, 4)


def compare_doscar_content(result, pdos=True):
    """Compare the content of a supplied get_quantity('poscar-structure') with reference data."""
    total_dos = np.array([(-3.44, -1.10400000e-43, -2.09900000e-43),
                          (-1.539, 1.40000000e-01, 2.66100000e-01),
                          (0.362, -3.62400000e-73, 2.00000000e+00),
                          (2.264, -1.33800000e-05, 2.00000000e+00),
                          (4.165, 3.15600000e+00, 8.00000000e+00),
                          (6.066, -2.41200000e-15, 8.00000000e+00),
                          (7.967, 3.15600000e+00, 1.40000000e+01),
                          (9.868, -1.38100000e-27, 1.40000000e+01),
                          (11.769, 2.90100000e+00, 1.95200000e+01),
                          (13.67, 0.00000000e+00, 2.00000000e+01)],
                         dtype=[('energy', '<f8'), ('total', '<f8'), ('integrated', '<f8')])  # yapf: disable
    partial_dos = np.array([[(-3.44 , -6.543e-45, -3.040e-46, -3.040e-46, -3.040e-46, 0., 0., 0., 0., 0.),
                             (-1.539,  8.295e-03,  3.853e-04,  3.853e-04,  3.853e-04, 0., 0., 0., 0., 0.),
                             ( 0.362, -1.336e-74, -3.400e-75, -3.400e-75, -3.400e-75, 0., 0., 0., 0., 0.),
                             ( 2.264, -4.932e-07, -1.255e-07, -1.255e-07, -1.255e-07, 0., 0., 0., 0., 0.),
                             ( 4.165,  1.163e-01,  2.961e-02,  2.961e-02,  2.961e-02, 0., 0., 0., 0., 0.),
                             ( 6.066, -5.870e-17, -2.973e-17, -2.973e-17, -2.973e-17, 0., 0., 0., 0., 0.),
                             ( 7.967,  7.681e-02,  3.890e-02,  3.890e-02,  3.890e-02, 0., 0., 0., 0., 0.),
                             ( 9.868, -1.988e-29, -1.936e-29, -1.923e-29, -1.944e-29, 0., 0., 0., 0., 0.),
                             (11.769,  5.486e-02,  3.117e-02,  3.938e-02,  3.827e-02, 0., 0., 0., 0., 0.),
                             (13.67 ,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00, 0., 0., 0., 0., 0.)],
                            [(-3.44 , -6.543e-45, -3.040e-46, -3.040e-46, -3.040e-46, 0., 0., 0., 0., 0.),
                             (-1.539,  8.295e-03,  3.853e-04,  3.853e-04,  3.853e-04, 0., 0., 0., 0., 0.),
                             ( 0.362, -1.336e-74, -3.400e-75, -3.400e-75, -3.400e-75, 0., 0., 0., 0., 0.),
                             ( 2.264, -4.932e-07, -1.255e-07, -1.255e-07, -1.255e-07, 0., 0., 0., 0., 0.),
                             ( 4.165,  1.163e-01,  2.961e-02,  2.961e-02,  2.961e-02, 0., 0., 0., 0., 0.),
                             ( 6.066, -5.870e-17, -2.973e-17, -2.973e-17, -2.973e-17, 0., 0., 0., 0., 0.),
                             ( 7.967,  7.681e-02,  3.890e-02,  3.890e-02,  3.890e-02, 0., 0., 0., 0., 0.),
                             ( 9.868, -1.988e-29, -1.936e-29, -1.923e-29, -1.944e-29, 0., 0., 0., 0., 0.),
                             (11.769,  5.486e-02,  3.117e-02,  3.938e-02,  3.827e-02, 0., 0., 0., 0., 0.),
                             (13.67 ,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00, 0., 0., 0., 0., 0.)],
                            [(-3.44 , -6.543e-45, -3.040e-46, -3.040e-46, -3.040e-46, 0., 0., 0., 0., 0.),
                          (-1.539,  8.295e-03,  3.853e-04,  3.853e-04,  3.853e-04, 0., 0., 0., 0., 0.),
                             ( 0.362, -1.336e-74, -3.400e-75, -3.400e-75, -3.400e-75, 0., 0., 0., 0., 0.),
                             ( 2.264, -4.932e-07, -1.255e-07, -1.255e-07, -1.255e-07, 0., 0., 0., 0., 0.),
                             ( 4.165,  1.163e-01,  2.961e-02,  2.961e-02,  2.961e-02, 0., 0., 0., 0., 0.),
                             ( 6.066, -5.870e-17, -2.973e-17, -2.973e-17, -2.973e-17, 0., 0., 0., 0., 0.),
                             ( 7.967,  7.681e-02,  3.890e-02,  3.890e-02,  3.890e-02, 0., 0., 0., 0., 0.),
                             ( 9.868, -1.988e-29, -1.936e-29, -1.923e-29, -1.944e-29, 0., 0., 0., 0., 0.),
                             (11.769,  5.486e-02,  3.117e-02,  3.938e-02,  3.827e-02, 0., 0., 0., 0., 0.),
                             (13.67 ,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00, 0., 0., 0., 0., 0.)],
                            [(-3.44 , -6.543e-45, -3.040e-46, -3.040e-46, -3.040e-46, 0., 0., 0., 0., 0.),
                             (-1.539,  8.295e-03,  3.853e-04,  3.853e-04,  3.853e-04, 0., 0., 0., 0., 0.),
                             ( 0.362, -1.336e-74, -3.400e-75, -3.400e-75, -3.400e-75, 0., 0., 0., 0., 0.),
                             ( 2.264, -4.932e-07, -1.255e-07, -1.255e-07, -1.255e-07, 0., 0., 0., 0., 0.),
                             ( 4.165,  1.163e-01,  2.961e-02,  2.961e-02,  2.961e-02, 0., 0., 0., 0., 0.),
                             ( 6.066, -5.870e-17, -2.973e-17, -2.973e-17, -2.973e-17, 0., 0., 0., 0., 0.),
                             ( 7.967,  7.681e-02,  3.890e-02,  3.890e-02,  3.890e-02, 0., 0., 0., 0., 0.),
                             ( 9.868, -1.988e-29, -1.936e-29, -1.923e-29, -1.944e-29, 0., 0., 0., 0., 0.),
                             (11.769,  5.486e-02,  3.117e-02,  3.938e-02,  3.827e-02, 0., 0., 0., 0., 0.),
                             (13.67 ,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00, 0., 0., 0., 0., 0.)]],
                           dtype=[('energy', '<f8'), ('s', '<f8'), ('py', '<f8'), ('px', '<f8'), ('pz', '<f8'), ('dxy', '<f8'), ('dyz', '<f8'), ('dz2', '<f8'), ('dxz', '<f8'), ('dx2-y2', '<f8')])
    result_total_dos = result['tdos']
    for item in total_dos.dtype.names:
        assert np.allclose(result_total_dos[item], total_dos[item])
    result_partial_dos = result['pdos']
    if pdos:
        result_partial_dos = result['pdos']
        for item in partial_dos.dtype.names:
            assert np.allclose(result_partial_dos[item], partial_dos[item])
    else:
        assert not result_partial_dos
    header = {'n_ions': 4, 'n_atoms': 4, 'cartesian': True, 'name': 'unknown system', 'emax': 13.67039808, 'emin': -3.43982524, 'n_dos': 10, 'efermi': 7.29482275, 'weight': 1.0}
    assert result['header'] == header
