"""Test the vasprun.xml parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy as np

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.aiida_utils import get_data_class


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_parse_vasprun(vasprun_parser):
    """Load a reference vasprun.xml and compare the result to a reference string."""

    quantity = vasprun_parser.get_quantity('occupancies')
    quantity = quantity['total']
    occ = quantity[0]
    occupancies = np.array([[[1., 1., 1., 1., 0.6667, 0.6667, 0.6667, -0., -0., -0.]]])
    assert occ.all() == occupancies.all()


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_parse_vasprun_version(vasprun_parser):
    """Load a reference vasprun.xml and fetch the VASP version."""
    version = vasprun_parser.get_quantity('version')
    assert version == '5.4.1'


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_parse_vasprun_kpoints(vasprun_parser):
    """Load a reference vasprun.xml and test that the parsed k-points are correct."""

    kpoints = vasprun_parser.get_quantity('kpoints')
    np.testing.assert_allclose(kpoints['points'][0], np.array([0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(kpoints['points'][-1], np.array([0.42857143, -0.42857143, 0.42857143]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_parse_vasprun_structure(vasprun_parser):
    """Load a reference vasprun.xml and test that the parsed k-points are correct."""
    structure = vasprun_parser.get_quantity('structure')
    # Check the unit cell
    np.testing.assert_allclose(structure['unitcell'][0], np.array([5.46503124, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(structure['unitcell'][1], np.array([0.0, 5.46503124, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(structure['unitcell'][2], np.array([0.0, 0.0, 5.46503124]), atol=0., rtol=1.0e-7)
    # Check first and last position
    np.testing.assert_allclose(structure['sites'][0]['position'], np.array([0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(structure['sites'][7]['position'], np.array([4.09877343, 4.09877343, 1.36625781]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_parse_vasprun_final_force(vasprun_parser):
    """Load a reference vasprun.xml and test that the forces are returned correctly."""
    forces = vasprun_parser.get_quantity('forces')
    forces = forces['final']
    forces_check = np.array([[-0.24286901, 0., 0.], [-0.24286901, 0., 0.], [3.41460162, 0., 0.], [0.44305748, 0., 0.],
                             [-0.73887169, 0.43727184, 0.43727184], [-0.94708885, -0.85011586, 0.85011586],
                             [-0.94708885, 0.85011586, -0.85011586], [-0.73887169, -0.43727184, -0.43727184]])
    # Check first, third and last position
    np.testing.assert_allclose(forces[0], forces_check[0], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(forces[2], forces_check[2], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(forces[7], forces_check[7], atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_parse_vasprun_final_stress(vasprun_parser):
    """Load a reference vasprun.xml and test that the stress are returned correctly."""
    stress = vasprun_parser.get_quantity('stress')
    stress = stress['final']
    stress_check = np.array([[-0.38703740, 0.00000000, 0.00000000], [0.00000000, 12.52362644, -25.93894358],
                             [0.00000000, -25.93894358, 12.52362644]])
    # Check entries
    np.testing.assert_allclose(stress[0], stress_check[0], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(stress[1], stress_check[1], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(stress[2], stress_check[2], atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('dielectric',)], indirect=True)
def test_parse_vasprun_dielectrics(vasprun_parser):
    """Load a reference vasprun.xml and test that the dielectrics are returned correctly."""
    dielectrics = vasprun_parser.get_quantity('dielectrics')
    imag = dielectrics['idiel']
    real = dielectrics['rdiel']
    energy = dielectrics['ediel']
    # Test shape of arrays
    assert imag.shape == (1000, 6)
    assert real.shape == (1000, 6)
    assert energy.shape == (1000,)
    # Test a few entries
    np.testing.assert_allclose(imag[0], np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(imag[500], np.array([0.0933, 0.0924, 0.0924, 0.0, 0.0082, 0.0]), atol=0., rtol=1.0e-7)

    np.testing.assert_allclose(imag[999], np.array([0.0035, 0.0035, 0.0035, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(real[0], np.array([12.0757, 11.4969, 11.4969, 0.0, 0.6477, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(real[500], np.array([-0.5237, -0.5366, -0.5366, 0.0, 0.0134, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(real[999],
                               np.array([6.57100000e-01, 6.55100000e-01, 6.55100000e-01, 0.0, -1.00000000e-04, 0.0]),
                               atol=0.,
                               rtol=1.0e-7)
    assert energy[500] == pytest.approx(10.2933)


@pytest.mark.parametrize(['vasprun_parser'], [('disp_details',)], indirect=True)
def test_parse_vasprun_epsilon(vasprun_parser):
    """Load a reference vasprun.xml and test that epsilon is returned correctly."""
    result = vasprun_parser.get_quantity('dielectrics')
    epsilon = result['epsilon']
    epsilon_ion = result['epsilon_ion']
    # Test shape of arrays
    assert epsilon.shape == (3, 3)
    assert epsilon_ion.shape == (3, 3)
    # Test a few entries
    test = np.array([[13.05544887, -0., 0.], [-0., 13.05544887, -0.], [0., 0., 13.05544887]])
    np.testing.assert_allclose(epsilon, test, atol=0., rtol=1.0e-7)
    test = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
    np.testing.assert_allclose(epsilon_ion, test, atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('localfield',)], indirect=True)
def test_parse_vasprun_born(vasprun_parser):
    """Load a reference vasprun.xml and test that the Born effective charges are
    returned correctly."""
    born = vasprun_parser.get_quantity('born_charges')
    born = born['born_charges']
    # Test shape of array
    assert born.shape == (8, 3, 3)
    # Test a few entries
    np.testing.assert_allclose(born[0][0], np.array([6.37225000e-03, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(born[0][-1], np.array([-4.21760000e-04, -2.19570210e-01, 3.20709600e-02]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(born[4][0], np.array([1.68565200e-01, -2.92058000e-02, -2.92058000e-02]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_parse_vasprun_dos(vasprun_parser):
    """Load a reference vasprun.xml and test that the density of states are
    returned correctly."""
    result = vasprun_parser.get_quantity('dos')
    dos = result['tdos']
    energy = result['energy']
    # Test shape of array
    assert dos.shape == (301,)
    assert energy.shape == (301,)
    # Test a few entries
    assert dos[150] == pytest.approx(4.1296)
    assert energy[150] == pytest.approx(2.3373)


@pytest.mark.parametrize(['vasprun_parser'], [('spin',)], indirect=True)
def test_parse_vasprun_dos_spin(vasprun_parser):
    """Load a reference vasprun.xml and test that the spin decomposed
    density of states is returned correctly."""
    result = vasprun_parser.get_quantity('dos')
    dos = result['tdos']
    # Test shape of array
    assert dos.shape == (
        2,
        1000,
    )
    # Test a few entries
    assert dos[0, 500] == pytest.approx(0.9839)
    assert dos[1, 500] == pytest.approx(0.9844)


@pytest.mark.parametrize(['vasprun_parser'], [('partial',)], indirect=True)
def test_parse_vasprun_pdos(vasprun_parser):
    """Load a reference vasprun.xml and test that the projected
    density of states is returned correctly."""
    result = vasprun_parser.get_quantity('dos')
    dos = result['pdos']
    energy = result['energy']
    # Test shape of array
    assert dos.shape == (8, 1000, 9)
    assert energy.shape == (1000,)
    # Test a few entries
    np.testing.assert_allclose(dos[3, 500], np.array([0.0770, 0.0146, 0.0109, 0.0155, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(dos[7, 500], np.array([0.0747, 0.0121, 0.0092, 0.0116, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    assert energy[500] == pytest.approx(0.01)


@pytest.mark.parametrize(['vasprun_parser'], [('partial',)], indirect=True)
def test_parse_vasprun_projectors(vasprun_parser):
    """Load a reference vasprun.xml and test that the state projectors are
    returned correctly."""
    proj = vasprun_parser.get_quantity('projectors')
    proj = proj['projectors']
    # Test shape of array
    assert proj.shape == (8, 64, 21, 9)
    # Test a few entries
    np.testing.assert_allclose(proj[0, 0, 5], np.array([0.0, 0.012, 0.0123, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(proj[7, 0, 5], np.array([0.1909, 0.0001, 0.0001, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(proj[4, 3, 5], np.array([0.2033, 0.0001, 0.0001, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_parse_vasprun_eigenvalues_occupancies(vasprun_parser):
    """Load a reference vasprun.xml and test that the eigenvalues are
    returned correctly."""
    eigen = vasprun_parser.get_quantity('eigenvalues')
    eigen = eigen['total']
    occ = vasprun_parser.get_quantity('occupancies')
    occ = occ['total']
    # Test shape of array
    assert eigen.shape == (64, 21)
    assert occ.shape == (64, 21)
    # Test a few entries
    assert eigen[0, 0] == pytest.approx(-6.2348)
    assert eigen[0, 15] == pytest.approx(5.8956)
    assert eigen[6, 4] == pytest.approx(-1.7424)
    assert occ[0, 0] == pytest.approx(1.0)
    assert occ[0, 15] == pytest.approx(0.6949)
    assert occ[6, 4] == pytest.approx(1.0)


@pytest.mark.parametrize(['vasprun_parser'], [('spin',)], indirect=True)
def test_parse_vasprun_eigenocc_spin_result(vasprun_parser):
    """Load a reference vasprun.xml and test that the spin decomposed eigenvalues
    are returned correctly."""
    eigen = vasprun_parser.get_quantity('eigenvalues')
    occ = vasprun_parser.get_quantity('occupancies')
    # Test shape of array
    assert eigen['up'].shape == (64, 25)
    assert occ['up'].shape == (64, 25)
    # Test a few entries
    assert eigen['up'][0, 0] == pytest.approx(-6.2363)
    assert eigen['up'][0, 15] == pytest.approx(5.8939)
    assert eigen['up'][6, 4] == pytest.approx(-1.7438)
    assert eigen['down'][0, 0] == pytest.approx(-6.2357)
    assert eigen['down'][0, 15] == pytest.approx(5.8946)
    assert eigen['down'][6, 4] == pytest.approx(-1.7432)
    assert occ['up'][0, 0] == pytest.approx(1.0)
    assert occ['up'][0, 15] == pytest.approx(0.6955)
    assert occ['up'][6, 4] == pytest.approx(1.0)
    assert occ['down'][0, 0] == pytest.approx(1.0)
    assert occ['down'][0, 15] == pytest.approx(0.6938)
    assert occ['down'][6, 4] == pytest.approx(1.0)


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_parse_vasprun_toten(vasprun_parser):
    """Load a reference vasprun.xml and test that one of the  total energies
    is returned correctly."""
    result = vasprun_parser.get_quantity('energies')
    assert set(result.keys()) == set(['energy_extrapolated_final', 'energy_extrapolated', 'electronic_steps'])
    energies = result['energy_extrapolated']
    test_array = np.array([-42.91113621])
    np.testing.assert_allclose(test_array, energies, atol=0., rtol=1.0e-7)
    # Test number of entries
    assert energies.shape == (1,)
    # Electronic steps should be one
    test_array = np.array([1])
    np.testing.assert_allclose(test_array, result['electronic_steps'], atol=0., rtol=1.0e-7)
    # Testing on VASP 5 so final total energy should not be the same as the last electronic step total energy.
    test_array = np.array([-0.00236711])
    np.testing.assert_allclose(test_array, result['energy_extrapolated_final'], atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [(['basic', 'vasprun.xml', {
    'energy_type': ['energy_free', 'energy_no_entropy']
}],)],
                         indirect=True)
def test_toten_multiple(vasprun_parser):
    """Load a reference vasprun.xml and test that multiple total energies
    are returned properly."""
    result = vasprun_parser.get_quantity('energies')
    assert set(result.keys()) == set(
        ['electronic_steps', 'energy_free', 'energy_free_final', 'energy_no_entropy', 'energy_no_entropy_final'])
    test_array = np.array([-42.91231976])
    np.testing.assert_allclose(test_array, result['energy_free'], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(test_array, result['energy_free_final'], atol=0., rtol=1.0e-7)
    test_array = np.array([-42.90995265])
    np.testing.assert_allclose(test_array, result['energy_no_entropy'], atol=0., rtol=1.0e-7)
    test_array = np.array([-42.91113621])
    np.testing.assert_allclose(test_array, result['energy_no_entropy_final'], atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [(['basic', 'vasprun.xml', {'electronic_step_energies': True}],)], indirect=True)
def test_parse_vasprun_toten_electronic(vasprun_parser):
    """Load a reference vasprun.xml and test that the total energies
    are returned correctly for the electronic steps."""
    result = vasprun_parser.get_quantity('energies')
    # Test that the default arrays are present
    assert set(result.keys()) == set(['energy_extrapolated_final', 'energy_extrapolated', 'electronic_steps'])
    energies = result['energy_extrapolated']
    test_array = np.array([-42.91113666, -42.91113621])
    np.testing.assert_allclose(test_array, energies, atol=0., rtol=1.0e-7)
    # Test number of entries
    assert energies.shape == (2,)
    # Electronic steps should be two
    test_array = np.array([2])
    np.testing.assert_allclose(test_array, result['electronic_steps'], atol=0., rtol=1.0e-7)
    # Testing on VASP 5 so final total energy should not be the same as the last electronic step total energy.
    test_array = np.array([-0.00236711])
    np.testing.assert_allclose(test_array, result['energy_extrapolated_final'], atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('relax',)], indirect=True)
def test_parse_vasprun_toten_relax(vasprun_parser):
    """Load a reference vasprun.xml and check that the total energies are
    returned correctly for relaxation runs."""
    result = vasprun_parser.get_quantity('energies')
    assert set(result.keys()) == set(['energy_extrapolated_final', 'energy_extrapolated', 'electronic_steps'])
    energies = result['energy_extrapolated']
    test_array = np.array([
        -42.91113348, -43.27757545, -43.36648855, -43.37734069, -43.38062479, -43.38334165, -43.38753003, -43.38708193, -43.38641449,
        -43.38701639, -43.38699488, -43.38773717, -43.38988315, -43.3898822, -43.39011239, -43.39020751, -43.39034244, -43.39044584,
        -43.39087657
    ])
    # Test energies
    np.testing.assert_allclose(test_array, energies, atol=0., rtol=1.0e-7)
    # Test number of entries
    assert energies.shape == test_array.shape
    # Electronic steps should be entries times one
    np.testing.assert_allclose(np.ones(19, dtype=int), result['electronic_steps'], atol=0., rtol=1.0e-7)
    # Testing on VASP 5 so final total energy should not be the same as the last electronic step total energy.
    test_array = np.array([
        -0.00236637, -0.00048614, -0.00047201, -0.00043261, -0.00041668, -0.00042584, -0.00043637, -0.00042806, -0.00042762, -0.00043875,
        -0.00042731, -0.00042705, -0.00043064, -0.00043051, -0.00043161, -0.00043078, -0.00043053, -0.00043149, -0.00043417
    ])
    np.testing.assert_allclose(test_array, result['energy_extrapolated_final'], atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [(['relax', 'vasprun.xml', {'electronic_step_energies': True}],)], indirect=True)
def test_parse_vasprun_toten_relax_electronic(vasprun_parser):
    """Load a reference vasprun.xml and check that the total energies
    are returned correctly for both the electronic and ionic steps."""
    result = vasprun_parser.get_quantity('energies')
    assert set(result.keys()) == set(['energy_extrapolated_final', 'energy_extrapolated', 'electronic_steps'])
    energies = result['energy_extrapolated']
    test_array_energies = [
        np.array([
            163.37398579, 14.26925896, -23.05190509, -34.91615104, -40.20080347, -42.18390876, -42.97469852, -43.31556073, -43.60169068,
            -43.61723125, -43.61871511, -43.61879751, -43.12548175, -42.90647187, -42.91031846, -42.91099027, -42.91111107, -42.91113348
        ]),
        np.array([-43.34236449, -43.31102002, -43.27768275, -43.27791002, -43.27761357, -43.27757545]),
        np.array([-43.40320524, -43.38084022, -43.36835045, -43.36666248, -43.36666583, -43.36649036, -43.36648855]),
        np.array([-43.37749056, -43.37749102, -43.37734414, -43.37734069]),
        np.array([-43.38117265, -43.38082881, -43.38063293, -43.38062479]),
        np.array([-43.38337336, -43.38334165]),
        np.array([-43.38778922, -43.38766017, -43.38752953, -43.38753003]),
        np.array([-43.38714489, -43.38708193]),
        np.array([-43.38640951, -43.38641449]),
        np.array([-43.3874799, -43.3871553, -43.38701949, -43.38701639]),
        np.array([-43.38790942, -43.38727062, -43.38700335, -43.38699488]),
        np.array([-43.38774394, -43.38773717]),
        np.array([-43.38984942, -43.3899134, -43.38988315]),
        np.array([-43.38988117, -43.3898822]),
        np.array([-43.39032165, -43.39017866, -43.39011239]),
        np.array([-43.39021044, -43.39020751]),
        np.array([-43.39034135, -43.39034244]),
        np.array([-43.39044466, -43.39044584]),
        np.array([-43.39084354, -43.39088709, -43.39087657])
    ]
    test_array_steps = np.array([18, 6, 7, 4, 4, 2, 4, 2, 2, 4, 4, 2, 3, 2, 3, 2, 2, 2, 3])
    # Build a flattened array (not using flatten from NumPy as the content is staggered) and
    # test number of electronic steps per ionic step
    test_array_energies_flattened = np.array([])
    for ionic_step in test_array_energies:
        test_array_energies_flattened = np.append(test_array_energies_flattened, ionic_step)
    assert energies.shape == test_array_energies_flattened.shape
    np.testing.assert_allclose(test_array_energies_flattened, energies, atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(test_array_steps, result['electronic_steps'], atol=0., rtol=1.0e-7)
    test_array_energies = np.array([
        -0.00236637, -0.00048614, -0.00047201, -0.00043261, -0.00041668, -0.00042584, -0.00043637, -0.00042806, -0.00042762, -0.00043875,
        -0.00042731, -0.00042705, -0.00043064, -0.00043051, -0.00043161, -0.00043078, -0.00043053, -0.00043149, -0.00043417
    ])
    # Testing on VASP 5 so final total energy should not be the same as the last electronic step total energy.
    np.testing.assert_allclose(test_array_energies, result['energy_extrapolated_final'], atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('disp',)], indirect=True)
def test_parse_vasprun_hessian(vasprun_parser):
    """Load a reference vasprun.xml and check that the Hessian matrix
    are returned correctly."""
    hessian = vasprun_parser.get_quantity('hessian')
    hessian = hessian['hessian']
    # Test shape
    assert hessian.shape == (24, 24)
    # Test a few entries
    assert np.allclose(
        hessian[0],
        np.array([
            -4.63550410e-01, 0.00000000e+00, 0.00000000e+00, -5.91774100e-02, 0.00000000e+00, 0.00000000e+00, 3.09711000e-02,
            0.00000000e+00, 0.00000000e+00, 3.20435400e-02, 0.00000000e+00, 0.00000000e+00, 1.15129840e-01, -8.16138200e-02, 8.17234700e-02,
            1.14879520e-01, 8.11324800e-02, 8.27409500e-02, 1.14879520e-01, -8.11324800e-02, -8.27409500e-02, 1.15129840e-01,
            8.16138200e-02, -8.17234700e-02
        ]))
    assert np.allclose(
        hessian[-2],
        np.array([
            8.16138200e-02, 1.15195590e-01, -8.38411100e-02, -8.17234700e-02, 1.14875090e-01, -8.53388100e-02, 3.46686900e-02,
            7.00672700e-02, 2.54288300e-02, -8.26222700e-02, 1.16185510e-01, 7.95575600e-02, -3.05970000e-04, 3.16827300e-02,
            2.86379000e-03, 5.42080000e-04, 3.27613500e-02, 1.12576000e-03, -1.34305000e-03, -5.86811100e-02, 2.83374000e-03,
            4.91688400e-02, -4.22101090e-01, 5.73736900e-02
        ]))


@pytest.mark.parametrize(['vasprun_parser'], [('disp',)], indirect=True)
def test_parse_vasprun_dynmat(vasprun_parser):
    """Load a reference vasprun.xml and check that the dynamical eigenvectors and eigenvalues
    are returned correctly."""
    result = vasprun_parser.get_quantity('dynmat')
    dynvec = result['dynvec']
    dyneig = result['dyneig']
    # test shape
    assert dynvec.shape == (24, 24)
    assert dyneig.shape == (24,)
    # test a few entries
    assert np.allclose(
        dynvec[0],
        np.array([
            7.28517310e-17, 7.25431601e-02, -4.51957676e-02, 1.15412776e-16, 4.51957676e-02, -7.25431601e-02, -1.37347223e-16,
            5.16257351e-01, -5.16257351e-01, 8.16789156e-17, 8.95098005e-02, -8.95098005e-02, -4.43838008e-17, -6.38031134e-02,
            6.38031134e-02, -1.80132830e-01, -2.97969516e-01, 2.97969516e-01, 1.80132830e-01, -2.97969516e-01, 2.97969516e-01,
            -2.09989969e-16, -6.38031134e-02, 6.38031134e-02
        ]))
    assert np.allclose(
        dynvec[4],
        np.array([
            -5.29825122e-13, -2.41759046e-01, -3.28913434e-01, -5.30734671e-13, -3.28913434e-01, -2.41759046e-01, 3.26325910e-13,
            -3.80807441e-02, -3.80807441e-02, -9.22956103e-13, -2.99868012e-01, -2.99868012e-01, 1.64418993e-01, 1.81002749e-01,
            1.81002749e-01, 3.11984195e-13, 2.73349550e-01, 2.73349550e-01, 2.59853610e-13, 2.73349550e-01, 2.73349550e-01, -1.64418993e-01,
            1.81002749e-01, 1.81002749e-01
        ]))
    assert dyneig[0] == pytest.approx(-1.36621537e+00)
    assert dyneig[4] == pytest.approx(-8.48939361e-01)
