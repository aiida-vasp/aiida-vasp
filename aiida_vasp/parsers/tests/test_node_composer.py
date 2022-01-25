"""Tests for NodeComposer."""
# pylint: disable=unused-wildcard-import,wildcard-import, redefined-outer-name, unused-argument
# pylint: disable=invalid-name
import pytest
import numpy as np

from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.utils.fixtures import *
from aiida_vasp.parsers.node_composer import NodeComposer, NODES_TYPES, clean_nan_values
from aiida_vasp.parsers.settings import NODES


@pytest.mark.parametrize(['poscar_parser'], [('poscar',)], indirect=True)
@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_create_node_structure(fresh_aiida_env, vasp_structure, poscar_parser):
    """Check that the node composer works for a StructureData node and contains the correct entries."""

    node_settings_key = 'poscar-structure'
    assert NODES[node_settings_key]['link_name'] == 'structure'
    assert NODES[node_settings_key]['type'] == 'structure'

    # Compose nodes
    composed_nodes = node_composer_test_helper('poscar-structure', NODES, [poscar_parser])

    # Compare
    result = composed_nodes.successful[NODES[node_settings_key]['link_name']]
    assert result.cell == vasp_structure.cell
    assert result.get_site_kindnames() == vasp_structure.get_site_kindnames()
    assert result.sites[2].position == vasp_structure.sites[2].position


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_create_node_kpoints(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for a KpointsData node and contains the correct entries."""

    node_settings_key = 'kpoints'
    assert NODES[node_settings_key]['link_name'] == 'kpoints'
    assert NODES[node_settings_key]['type'] == 'array.kpoints'

    # Compose nodes
    composed_nodes = node_composer_test_helper('kpoints', NODES, [vasprun_parser])

    # Compare
    kpoints = composed_nodes.successful['kpoints']
    np.testing.assert_allclose(kpoints.get_kpoints()[0], np.array([0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(kpoints.get_kpoints()[-1], np.array([0.42857143, -0.42857143, 0.42857143]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['outcar_parser'], [(['disp_details', 'OUTCAR'],)], indirect=True)
def test_create_node_misc(fresh_aiida_env, outcar_parser):
    """Check that the node composer works for the default misc node and contains the correct entries."""

    # Prepare inputs for the node composer
    quantities_to_parse = ['run_stats', 'run_status']
    requested_node = {'misc': {'link_name': 'misc', 'type': 'dict', 'quantities': quantities_to_parse}}
    # Compose nodes
    composed_nodes = node_composer_test_helper('misc', requested_node, [outcar_parser])

    # Compare
    data_dict = composed_nodes.successful['misc'].get_dict()
    assert data_dict['run_stats']
    assert data_dict['run_stats']['total_cpu_time_used'] == pytest.approx(89.795)
    assert data_dict['run_stats']['average_memory_used'] == pytest.approx(0.0)

    assert data_dict['run_status']['last_iteration_index'] == [15, 5]
    assert data_dict['run_status']['finished']
    assert data_dict['run_status']['ionic_converged']
    assert data_dict['run_status']['electronic_converged']
    assert data_dict['run_status']['nelm'] == 60
    assert data_dict['run_status']['nsw'] == 61


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_create_node_forces(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the forces node and contains the correct forces."""

    node_settings_key = 'forces'
    assert NODES[node_settings_key]['link_name'] == 'forces'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    forces_check = np.array([[-0.24286901, 0., 0.], [-0.24286901, 0., 0.], [3.41460162, 0., 0.], [0.44305748, 0., 0.],
                             [-0.73887169, 0.43727184, 0.43727184], [-0.94708885, -0.85011586, 0.85011586],
                             [-0.94708885, 0.85011586, -0.85011586], [-0.73887169, -0.43727184, -0.43727184]])
    data_obj = composed_nodes.successful['forces']
    forces = data_obj.get_array('final')
    # First, third and last position
    np.testing.assert_allclose(forces[0], forces_check[0], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(forces[2], forces_check[2], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(forces[7], forces_check[7], atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_create_node_stress(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the stress node and contains the correct stress."""

    node_settings_key = 'stress'
    assert NODES[node_settings_key]['link_name'] == 'stress'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    stress_check = np.array([[-0.38703740, 0.00000000, 0.00000000], [0.00000000, 12.52362644, -25.93894358],
                             [0.00000000, -25.93894358, 12.52362644]])
    data_obj = composed_nodes.successful['stress']
    stress = data_obj.get_array('final')
    # First, third and last position
    np.testing.assert_allclose(stress[0], stress_check[0], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(stress[1], stress_check[1], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(stress[2], stress_check[2], atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_create_node_trajectory_forces(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the trajectory node and contain the correct forces."""

    node_settings_key = 'trajectory'
    assert NODES[node_settings_key]['link_name'] == 'trajectory'
    assert NODES[node_settings_key]['type'] == 'array.trajectory'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['trajectory']
    traj = data_obj.get_array('forces')
    np.testing.assert_allclose(traj[0][0], np.array([-0.24286901, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(traj[0][-1], np.array([-0.73887169, -0.43727184, -0.43727184]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(traj[0][-1], traj[1][-1], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(traj[0][0], traj[1][0], atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('relax',)], indirect=True)
def test_create_node_trajectory_forces_ionic(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the trajectory node and contain the correct forces when ionic steps is present."""

    node_settings_key = 'trajectory'
    assert NODES[node_settings_key]['link_name'] == 'trajectory'
    assert NODES[node_settings_key]['type'] == 'array.trajectory'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['trajectory']
    traj = data_obj.get_array('forces')
    assert traj.shape == (19, 8, 3)
    # First and last atom
    np.testing.assert_allclose(traj[0][0], np.array([-2.42632080e-01, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(traj[0][-1], np.array([-7.38879520e-01, -4.37063010e-01, -4.37063010e-01]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(traj[-1][0], np.array([1.55852000e-03, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(traj[-1][-1], np.array([-1.75970000e-03, 1.12150000e-04, 1.12150000e-04]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('relax',)], indirect=True)
def test_create_node_trajectory_cells(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the trajectory node and contain the correct cells when ionic steps is present."""

    node_settings_key = 'trajectory'
    assert NODES[node_settings_key]['link_name'] == 'trajectory'
    assert NODES[node_settings_key]['type'] == 'array.trajectory'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['trajectory']
    traj = data_obj.get_array('cells')
    assert traj.shape == (19, 3, 3)
    # First and last entry
    np.testing.assert_allclose(traj[0][0], np.array([5.46503124e+00, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(traj[0][-1], np.array([0.0, 0.0, 5.46503124e+00]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(traj[-1][0], np.array([5.46702248e+00, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(traj[-1][-1], np.array([0.0, 2.19104000e-03, 5.46705225e+00]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('relax',)], indirect=True)
def test_create_node_trajectory_positions(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the trajectory node and contain the correct positions when ionic steps is present."""

    node_settings_key = 'trajectory'
    assert NODES[node_settings_key]['link_name'] == 'trajectory'
    assert NODES[node_settings_key]['type'] == 'array.trajectory'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['trajectory']
    traj = data_obj.get_array('positions')
    assert traj.shape == (19, 8, 3)
    # First and last entry
    np.testing.assert_allclose(traj[0][0], np.array([0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(traj[0][-1], np.array([0.75, 0.75, 0.25]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(traj[-1][0], np.array([-0.00621692, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(traj[-1][-1], np.array([0.7437189, 0.74989833, 0.24989833]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('dielectric',)], indirect=True)
def test_create_node_dielectrics(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the dielectrics node and contain the correct dielectrics."""

    node_settings_key = 'dielectrics'
    assert NODES[node_settings_key]['link_name'] == 'dielectrics'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['dielectrics']
    imag = data_obj.get_array('idiel')
    real = data_obj.get_array('rdiel')
    energy = data_obj.get_array('ediel')
    # Shape of arrays
    assert imag.shape == (1000, 6)
    assert real.shape == (1000, 6)
    assert energy.shape == (1000,)
    # A few entries
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
def test_create_node_dielectrics_epsilon(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the dielectrics node and contain the correct epsilon."""

    node_settings_key = 'dielectrics'
    assert NODES[node_settings_key]['link_name'] == 'dielectrics'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['dielectrics']
    epsilon = data_obj.get_array('epsilon')
    epsilon_ion = data_obj.get_array('epsilon_ion')
    # Shape of arrays
    assert epsilon.shape == (3, 3)
    assert epsilon_ion.shape == (3, 3)
    # A few entries
    test = np.array([[13.05544887, -0., 0.], [-0., 13.05544887, -0.], [0., 0., 13.05544887]])
    np.testing.assert_allclose(epsilon, test, atol=0., rtol=1.0e-7)
    test = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
    np.testing.assert_allclose(epsilon_ion, test, atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('localfield',)], indirect=True)
def test_create_node_born_charges(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the Born charges node and contain the correct Born charges."""

    node_settings_key = 'born_charges'
    assert NODES[node_settings_key]['link_name'] == 'born_charges'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['born_charges']
    born = data_obj.get_array('born_charges')
    # Shape of array
    assert born.shape == (8, 3, 3)
    # A few entries
    np.testing.assert_allclose(born[0][0], np.array([6.37225000e-03, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(born[0][-1], np.array([-4.21760000e-04, -2.19570210e-01, 3.20709600e-02]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(born[4][0], np.array([1.68565200e-01, -2.92058000e-02, -2.92058000e-02]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_create_node_dos(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the density of states node and contain the correct density of states."""

    node_settings_key = 'dos'
    assert NODES[node_settings_key]['link_name'] == 'dos'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['dos']
    dos = data_obj.get_array('tdos')
    energy = data_obj.get_array('energy')
    # Shape of array
    assert dos.shape == (301,)
    assert energy.shape == (301,)
    # A few entries
    assert dos[150] == pytest.approx(4.1296)
    assert energy[150] == pytest.approx(2.3373)


@pytest.mark.parametrize(['vasprun_parser'], [('spin',)], indirect=True)
def test_create_node_dos_spin(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the density of states node and contain the correct density of states for spins."""

    node_settings_key = 'dos'
    assert NODES[node_settings_key]['link_name'] == 'dos'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['dos']
    dos = data_obj.get_array('tdos')
    # Shape of array
    assert dos.shape == (
        2,
        1000,
    )
    # A few entries
    assert dos[0, 500] == pytest.approx(0.9839)
    assert dos[1, 500] == pytest.approx(0.9844)


@pytest.mark.parametrize(['vasprun_parser'], [('partial',)], indirect=True)
def test_create_node_dos_partial(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the density of states node and contain the correct decomposed density of states."""

    node_settings_key = 'dos'
    assert NODES[node_settings_key]['link_name'] == 'dos'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['dos']
    dos = data_obj.get_array('pdos')
    energy = data_obj.get_array('energy')
    # Shape of array
    assert dos.shape == (8, 1000, 9)
    assert energy.shape == (1000,)
    # A few entries
    np.testing.assert_allclose(dos[3, 500], np.array([0.0770, 0.0146, 0.0109, 0.0155, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(dos[7, 500], np.array([0.0747, 0.0121, 0.0092, 0.0116, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    assert energy[500] == pytest.approx(0.01)


@pytest.mark.parametrize(['vasprun_parser'], [('partial',)], indirect=True)
def test_create_node_projectors(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the projectors node and contain the correct projectors."""

    node_settings_key = 'projectors'
    assert NODES[node_settings_key]['link_name'] == 'projectors'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['projectors']
    proj = data_obj.get_array('projectors')
    # Shape of array
    assert proj.shape == (8, 64, 21, 9)
    # A few entries
    np.testing.assert_allclose(proj[0, 0, 5], np.array([0.0, 0.012, 0.0123, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(proj[7, 0, 5], np.array([0.1909, 0.0001, 0.0001, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(proj[4, 3, 5], np.array([0.2033, 0.0001, 0.0001, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_create_node_bands(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the bands node and contain the correct eigenvalues and occupancies."""

    node_settings_key = 'bands'
    assert NODES[node_settings_key]['link_name'] == 'bands'
    assert NODES[node_settings_key]['type'] == 'array.bands'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['bands']
    eigenocc = data_obj.get_bands(also_occupations=True)
    eigen = eigenocc[0]
    occ = eigenocc[1]
    # Shape of array
    assert eigen.shape == (1, 64, 21)
    assert occ.shape == (1, 64, 21)
    # A few entries
    assert eigen[0, 0, 0] == pytest.approx(-6.2348)
    assert eigen[0, 0, 15] == pytest.approx(5.8956)
    assert eigen[0, 6, 4] == pytest.approx(-1.7424)
    assert occ[0, 0, 0] == pytest.approx(1.0)
    assert occ[0, 0, 15] == pytest.approx(0.6949)
    assert occ[0, 6, 4] == pytest.approx(1.0)


@pytest.mark.parametrize(['vasprun_parser'], [('spin',)], indirect=True)
def test_create_node_bands_spin(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the bands node and contain the correct eigenvalues and occupancies for spin."""

    node_settings_key = 'bands'
    assert NODES[node_settings_key]['link_name'] == 'bands'
    assert NODES[node_settings_key]['type'] == 'array.bands'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['bands']
    eigenocc = data_obj.get_bands(also_occupations=True)
    eigen = eigenocc[0]
    occ = eigenocc[1]
    # Shape of array
    assert eigen.shape == (2, 64, 25)
    assert occ.shape == (2, 64, 25)
    # A few entries
    assert eigen[0, 0, 0] == pytest.approx(-6.2363)
    assert eigen[0, 0, 15] == pytest.approx(5.8939)
    assert eigen[0, 6, 4] == pytest.approx(-1.7438)
    assert eigen[1, 0, 0] == pytest.approx(-6.2357)
    assert eigen[1, 0, 15] == pytest.approx(5.8946)
    assert eigen[1, 6, 4] == pytest.approx(-1.7432)
    assert occ[0, 0, 0] == pytest.approx(1.0)
    assert occ[0, 0, 15] == pytest.approx(0.6955)
    assert occ[0, 6, 4] == pytest.approx(1.0)
    assert occ[1, 0, 0] == pytest.approx(1.0)
    assert occ[1, 0, 15] == pytest.approx(0.6938)
    assert occ[1, 6, 4] == pytest.approx(1.0)


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_create_node_energies(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the energies node and contain the correct energies."""

    node_settings_key = 'energies'
    assert NODES[node_settings_key]['link_name'] == 'energies'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    energies = composed_nodes.successful['energies']
    assert set(energies.get_arraynames()) == set(['energy_extrapolated_final', 'energy_extrapolated', 'electronic_steps'])
    energies_ext = energies.get_array('energy_extrapolated')
    test_array = np.array([-42.91113621])
    np.testing.assert_allclose(test_array, energies_ext, atol=0., rtol=1.0e-7)
    # Test number of entries
    assert energies_ext.shape == (1,)
    # Electronic steps should be one
    test_array = np.array([1])
    np.testing.assert_allclose(test_array, energies.get_array('electronic_steps'), atol=0., rtol=1.0e-7)
    # Testing on VASP 5 so final total energy should not be the same as the last electronic step total energy.
    test_array = np.array([-0.00236711])
    np.testing.assert_allclose(test_array, energies.get_array('energy_extrapolated_final'), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [(['basic', 'vasprun.xml', {
    'energy_type': ['energy_free', 'energy_no_entropy']
}],)],
                         indirect=True)
def test_create_node_energies_multiple(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the energies node and contain the correct set of chosen energies."""

    node_settings_key = 'energies'
    assert NODES[node_settings_key]['link_name'] == 'energies'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    energies = composed_nodes.successful['energies']
    assert set(energies.get_arraynames()) == set(
        ['electronic_steps', 'energy_free', 'energy_free_final', 'energy_no_entropy', 'energy_no_entropy_final'])
    test_array = np.array([-42.91231976])
    np.testing.assert_allclose(test_array, energies.get_array('energy_free'), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(test_array, energies.get_array('energy_free_final'), atol=0., rtol=1.0e-7)
    test_array = np.array([-42.90995265])
    np.testing.assert_allclose(test_array, energies.get_array('energy_no_entropy'), atol=0., rtol=1.0e-7)
    test_array = np.array([-42.91113621])
    np.testing.assert_allclose(test_array, energies.get_array('energy_no_entropy_final'), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [(['basic', 'vasprun.xml', {'electronic_step_energies': True}],)], indirect=True)
def test_create_node_energies_electronic(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the energies node and contain the correct energies for each electronic step update"""

    node_settings_key = 'energies'
    assert NODES[node_settings_key]['link_name'] == 'energies'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    energies = composed_nodes.successful['energies']
    energies_ext = energies.get_array('energy_extrapolated')
    test_array = np.array([-42.91113666, -42.91113621])
    np.testing.assert_allclose(test_array, energies_ext, atol=0., rtol=1.0e-7)
    # Test number of entries
    assert energies_ext.shape == (2,)
    # Electronic steps should be two
    test_array = np.array([2])
    np.testing.assert_allclose(test_array, energies.get_array('electronic_steps'), atol=0., rtol=1.0e-7)
    # Testing on VASP 5 so final total energy should not be the same as the last electronic step total energy.
    test_array = np.array([-0.00236711])
    np.testing.assert_allclose(test_array, energies.get_array('energy_extrapolated_final'), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('relax',)], indirect=True)
def test_create_node_energies_relax(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the energies node and contain the correct energies for each ionic step."""

    node_settings_key = 'energies'
    assert NODES[node_settings_key]['link_name'] == 'energies'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    energies = composed_nodes.successful['energies']
    energies_ext = energies.get_array('energy_extrapolated')
    assert set(energies.get_arraynames()) == set(['energy_extrapolated_final', 'energy_extrapolated', 'electronic_steps'])
    test_array = np.array([
        -42.91113348, -43.27757545, -43.36648855, -43.37734069, -43.38062479, -43.38334165, -43.38753003, -43.38708193, -43.38641449,
        -43.38701639, -43.38699488, -43.38773717, -43.38988315, -43.3898822, -43.39011239, -43.39020751, -43.39034244, -43.39044584,
        -43.39087657
    ])
    # Test energies
    np.testing.assert_allclose(test_array, energies_ext, atol=0., rtol=1.0e-7)
    # Test number of entries
    assert energies_ext.shape == test_array.shape
    # Electronic steps should be entries times one
    np.testing.assert_allclose(np.ones(19, dtype=int), energies.get_array('electronic_steps'), atol=0., rtol=1.0e-7)
    # Testing on VASP 5 so final total energy should not be the same as the last electronic step total energy.
    test_array = np.array([
        -0.00236637, -0.00048614, -0.00047201, -0.00043261, -0.00041668, -0.00042584, -0.00043637, -0.00042806, -0.00042762, -0.00043875,
        -0.00042731, -0.00042705, -0.00043064, -0.00043051, -0.00043161, -0.00043078, -0.00043053, -0.00043149, -0.00043417
    ])
    np.testing.assert_allclose(test_array, energies.get_array('energy_extrapolated_final'), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [(['relax', 'vasprun.xml', {'electronic_step_energies': True}],)], indirect=True)
def test_create_node_energies_electronic_relax(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the energies node and contain the correct energies for each ionic and ionic step."""

    node_settings_key = 'energies'
    assert NODES[node_settings_key]['link_name'] == 'energies'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    energies = composed_nodes.successful['energies']
    energies_ext = energies.get_array('energy_extrapolated')
    assert set(energies.get_arraynames()) == set(['energy_extrapolated_final', 'energy_extrapolated', 'electronic_steps'])
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
    assert energies_ext.shape == test_array_energies_flattened.shape
    np.testing.assert_allclose(test_array_energies_flattened, energies_ext, atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(test_array_steps, energies.get_array('electronic_steps'), atol=0., rtol=1.0e-7)
    test_array_energies = np.array([
        -0.00236637, -0.00048614, -0.00047201, -0.00043261, -0.00041668, -0.00042584, -0.00043637, -0.00042806, -0.00042762, -0.00043875,
        -0.00042731, -0.00042705, -0.00043064, -0.00043051, -0.00043161, -0.00043078, -0.00043053, -0.00043149, -0.00043417
    ])
    # Testing on VASP 5 so final total energy should not be the same as the last electronic step total energy.
    np.testing.assert_allclose(test_array_energies, energies.get_array('energy_extrapolated_final'), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize(['vasprun_parser'], [('disp',)], indirect=True)
def test_create_node_hessian(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the hessian node and contain the correct hessian."""

    node_settings_key = 'hessian'
    assert NODES[node_settings_key]['link_name'] == 'hessian'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['hessian']
    hessian = data_obj.get_array('hessian')
    # Test shape
    assert hessian.shape == (24, 24)
    # A few entries
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
def test_create_node_dynmat(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for the dynmat node and contain the correct dynamical matrix."""

    node_settings_key = 'dynmat'
    assert NODES[node_settings_key]['link_name'] == 'dynmat'
    assert NODES[node_settings_key]['type'] == 'array'

    # Compose nodes
    composed_nodes = node_composer_test_helper(node_settings_key, NODES, [vasprun_parser])

    # Compare
    data_obj = composed_nodes.successful['dynmat']
    dynvec = data_obj.get_array('dynvec')
    dyneig = data_obj.get_array('dyneig')
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


@pytest.mark.parametrize(['outcar_parser'],
                         [(['disp_details', 'OUTCAR', {
                             'quantities_to_parse': ['symmetries', 'elastic_moduli', 'run_status', 'run_stats']
                         }],)],
                         indirect=True)
def test_create_node_dict_custom(fresh_aiida_env, outcar_parser):
    """Check that the node composer works for custom Dict node."""

    # Prepare inputs for the node composer
    quantities_to_parse = ['symmetries', 'elastic_moduli', 'run_stats', 'run_status']
    requested_node = {'misc': {'link_name': 'my_custom_node', 'type': 'dict', 'quantities': quantities_to_parse}}
    parsed_quantities = {}
    equivalent_keys = {}
    for item in quantities_to_parse:
        parsed_quantities[item] = outcar_parser.get_quantity(item)
        equivalent_keys[item] = [item]

    # Compose node
    composed_nodes = NodeComposer(requested_node, equivalent_keys, parsed_quantities)
    assert isinstance(composed_nodes.successful['my_custom_node'], get_data_class('dict'))

    # Compare
    compare_symmetries = {
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
    data_dict = composed_nodes.successful['my_custom_node'].get_dict()
    assert set(data_dict['symmetries']) == set(compare_symmetries)

    # then elastic moduli
    test = np.array([1674.5786, 704.739, 704.739, -0.0, 0.0, 0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['symmetrized'][0], test)
    test = np.array([0.0, 0.0, 0.0, -0.0, -0.0, 1122.6622])
    np.testing.assert_allclose(data_dict['elastic_moduli']['symmetrized'][5], test)
    test = np.array([705.0238, 1674.8491, 705.0238, -0.0, -0.0, 0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['non_symmetrized'][1], test)
    test = np.array([-0.0078, -0.0495, 0.0147, 0.0, 1123.0829, -0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['non_symmetrized'][4], test)
    test = np.array([704.739, 704.739, 1674.5786, -0.0, -0.0, 0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['total'][2], test)
    test = np.array([-0.0, -0.0, -0.0, 775.8054, 0.0, -0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['total'][3], test)

    assert data_dict['run_stats']
    assert data_dict['run_stats']['total_cpu_time_used'] == pytest.approx(89.795)
    assert data_dict['run_stats']['average_memory_used'] == pytest.approx(0.0)

    assert data_dict['run_status']['last_iteration_index'] == [15, 5]
    assert data_dict['run_status']['finished']
    assert data_dict['run_status']['ionic_converged']
    assert data_dict['run_status']['electronic_converged']
    assert data_dict['run_status']['nelm'] == 60
    assert data_dict['run_status']['nsw'] == 61


@pytest.mark.parametrize(
    ['vasprun_parser'],
    [(['basic', 'vasprun.xml', {
        'quantities_to_parse': ['fermi_level', 'total_energies', 'energies', 'maximum_force', 'maximum_stress']
    }],)],
    indirect=True)
def test_create_node_dict_custom2(fresh_aiida_env, vasprun_parser):
    """Check that the node composer works for yet another custom Dict node."""

    # Prepare inputs for the node composer
    quantities_to_parse = ['fermi_level', 'total_energies', 'energies', 'maximum_force', 'maximum_stress']
    requested_node = {'misc': {'link_name': 'my_custom_node', 'type': 'dict', 'quantities': quantities_to_parse}}
    parsed_quantities = {}
    equivalent_keys = {}
    for item in quantities_to_parse:
        parsed_quantities[item] = vasprun_parser.get_quantity(item)
        equivalent_keys[item] = [item]

    # Compose node
    composed_nodes = NodeComposer(requested_node, equivalent_keys, parsed_quantities)
    assert isinstance(composed_nodes.successful['my_custom_node'], get_data_class('dict'))

    # Compare
    data_dict = composed_nodes.successful['my_custom_node'].get_dict()
    assert data_dict['fermi_level'] == pytest.approx(5.96764939)
    assert data_dict['total_energies']['energy_extrapolated'] == pytest.approx(-42.91113621)
    assert data_dict['energies']['energy_extrapolated'][0] == pytest.approx(-42.91113621)
    assert data_dict['maximum_stress'] == pytest.approx(28.803993008871014)
    assert data_dict['maximum_force'] == pytest.approx(3.41460162)


@pytest.mark.parametrize(['outcar_parser'], [(['basic_run', 'OUTCAR', {'quantities_to_parse': NODES_TYPES['dict']}],)], indirect=True)
@pytest.mark.parametrize(['vasprun_parser'], [(['basic_run', 'vasprun.xml', {'quantities_to_parse': NODES_TYPES['dict']}],)], indirect=True)
def test_create_node_misc_all(fresh_aiida_env, vasprun_parser, outcar_parser):
    """Check that the node composer works for the misc node containing all valid quantities."""

    # Prepare inputs for the node composer
    quantities_to_parse = NODES_TYPES['dict']
    requested_node = {'misc': {'link_name': 'misc', 'type': 'dict', 'quantities': quantities_to_parse}}

    # Compose nodes
    composed_nodes = node_composer_test_helper('misc', requested_node, [vasprun_parser, outcar_parser])

    # Compare
    misc = composed_nodes.successful['misc'].get_dict()
    assert misc['band_properties']['cbm'] == pytest.approx(5.075)
    assert misc['band_properties']['vbm'] == pytest.approx(4.2811)
    assert misc['band_properties']['band_gap'] == pytest.approx(0.793899999999999)
    assert not misc['band_properties']['is_direct_gap']
    assert misc['version'] == '5.3.5'
    assert misc['total_energies']['energy_extrapolated'] == pytest.approx(-36.09616894)
    assert misc['maximum_stress'] == pytest.approx(8.50955439)
    assert misc['maximum_force'] == pytest.approx(0.0)
    assert misc['run_status']['nelm'] == 60
    assert misc['run_status']['last_iteration_index'] == [1, 2]
    assert misc['run_status']['nsw'] == 0
    assert misc['run_status']['finished']
    assert misc['run_status']['ionic_converged'] is None
    assert misc['run_status']['electronic_converged']
    assert not misc['run_status']['consistent_nelm_breach']
    assert not misc['run_status']['contains_nelm_breach']
    assert misc['run_stats']['mem_usage_base'] == pytest.approx(30000.0)
    assert misc['run_stats']['mem_usage_nonl-proj'] == pytest.approx(7219.0)
    assert misc['run_stats']['mem_usage_fftplans'] == pytest.approx(776.0)
    assert misc['run_stats']['mem_usage_grid'] == pytest.approx(1605.0)
    assert misc['run_stats']['mem_usage_one-center'] == pytest.approx(124.0)
    assert misc['run_stats']['mem_usage_wavefun'] == pytest.approx(814.0)
    assert misc['run_stats']['maximum_memory_used'] == pytest.approx(95344.0)
    assert misc['run_stats']['average_memory_used'] == pytest.approx(0.0)
    assert misc['run_stats']['total_cpu_time_used'] == pytest.approx(20.463)
    assert misc['run_stats']['user_time'] == pytest.approx(11.266)
    assert misc['run_stats']['system_time'] == pytest.approx(9.197)
    assert misc['run_stats']['elapsed_time'] == pytest.approx(22.518)
    assert 'symmetries' in misc
    # No magnetization
    assert misc['magnetization'] is None
    assert misc['site_magnetization'] == {
        'sphere': {
            'x': {
                'site_moment': {},
                'total_magnetization': {}
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
        'full_cell': []
    }


def test_nan_inf_cleaning():
    """
    Test cleaning nan/inf values
    """

    example = {'a': 1.0, 'b': {'c': 2.0, 'd': np.inf, 'e': np.nan}, 'd': {'e': {'g': np.inf}}}
    clean_nan_values(example)
    assert example['a'] == pytest.approx(1.0)
    assert example['b']['d'] == 'inf'
    assert example['b']['e'] == 'nan'
    assert example['d']['e']['g'] == 'inf'


def node_composer_test_helper(node_settings_key, NODES, parsers):
    """A helper routine to minimize code ducplication across tests for node composition."""

    # Prepare inputs for tqhe node composer
    requested_node = {NODES[node_settings_key]['link_name']: NODES[node_settings_key]}
    parsed_quantities = {}
    equivalent_keys = {}
    for parser in parsers:
        # Assume there are no similar quantities defined among different parsers. They should be prefixed
        # if not default.
        for item in NODES[node_settings_key]['quantities']:
            if item in parser.PARSABLE_QUANTITIES:
                parsed_quantities[item] = parser.get_quantity(item)
                equivalent_keys[item] = item

    # Compose node
    composed_nodes = NodeComposer(requested_node, equivalent_keys, parsed_quantities)
    data_class = get_data_class(NODES[node_settings_key]['type'])
    assert NODES[node_settings_key]['link_name'] in composed_nodes.successful
    assert isinstance(composed_nodes.successful[NODES[node_settings_key]['link_name']], data_class)

    return composed_nodes
