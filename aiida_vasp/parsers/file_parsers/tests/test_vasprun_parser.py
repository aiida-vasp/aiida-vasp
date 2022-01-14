"""Test the vasprun.xml parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy as np

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.parsers.node_composer import NodeComposer, get_node_composer_inputs_from_file_parser, clean_nan_values


def test_version(fresh_aiida_env, vasprun_parser):
    """Parse a reference vasprun.xml and fetch the VASP version."""
    quantity = vasprun_parser.get_quantity('version')
    assert quantity == '5.4.4'


def test_parse_vasprun(fresh_aiida_env, vasprun_parser):
    """Parse a reference vasprun.xml file with the VasprunParser and compare the result to a reference string."""

    quantity = vasprun_parser.get_quantity('occupancies')
    occ = quantity[0]
    occupancies = np.array([[[1., 1., 1., 1., 0.6667, 0.6667, 0.6667, -0., -0., -0.]]])
    assert occ.all() == occupancies.all()
    # eFL: How do we want to store scalar values?
    #assert  == 7.29482275


@pytest.mark.parametrize('vasprun_parser', [('basic', {})], indirect=True)
def test_parameter_results(fresh_aiida_env, vasprun_parser):
    """
    Test that the parameter node is a ParametersData instance.

    Should contain the Fermi level, total_energies, maximum_force and maximum_stress.

    """

    vasprun_parser._settings._output_nodes_dict.update({  # pylint: disable=protected-access
        'misc': {
            'type': 'dict',
            'quantities': ['fermi_level', 'total_energies', 'energies', 'maximum_force', 'maximum_stress'],
            'link_name': 'my_custom_node'
        }
    })

    quantity_keys = ['fermi_level', 'total_energies', 'energies', 'maximum_force', 'maximum_stress']
    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=quantity_keys)
    data_obj = NodeComposer.compose('dict', inputs)

    ref_class = get_data_class('dict')
    assert isinstance(data_obj, ref_class)
    data_dict = data_obj.get_dict()
    assert data_dict['fermi_level'] == pytest.approx(5.96764939)
    assert data_dict['total_energies']['energy_extrapolated'] == pytest.approx(-42.91113621)
    assert data_dict['energies']['energy_extrapolated'][0] == pytest.approx(-42.91113621)
    assert data_dict['maximum_stress'] == pytest.approx(28.803993008871014)
    assert data_dict['maximum_force'] == pytest.approx(3.41460162)


@pytest.mark.parametrize('vasprun_parser', [('basic', {})], indirect=True)
def test_kpoints(fresh_aiida_env, vasprun_parser):
    """Test that the kpoints result node is a KpointsData instance."""

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['kpoints'])
    data_obj = NodeComposer.compose('array.kpoints', inputs)

    ref_class = get_data_class('array.kpoints')
    assert isinstance(data_obj, ref_class)
    np.testing.assert_allclose(data_obj.get_kpoints()[0], np.array([0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj.get_kpoints()[-1], np.array([0.42857143, -0.42857143, 0.42857143]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('basic', {})], indirect=True)
def test_structure(fresh_aiida_env, vasprun_parser):
    """
    Test that the structure result node is a StructureData instance.

    Also check various other important properties.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['structure'])
    data_obj = NodeComposer.compose('structure', inputs)
    # check object
    ref_obj = get_data_class('structure')
    assert isinstance(data_obj, ref_obj)
    # check cell
    unitcell = data_obj.cell
    np.testing.assert_allclose(unitcell[0], np.array([5.46503124, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(unitcell[1], np.array([0.0, 5.46503124, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(unitcell[2], np.array([0.0, 0.0, 5.46503124]), atol=0., rtol=1.0e-7)
    # check first and last position
    np.testing.assert_allclose(data_obj.sites[0].position, np.array([0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj.sites[7].position, np.array([4.09877343, 4.09877343, 1.36625781]), atol=0., rtol=1.0e-7)
    # check volume
    assert data_obj.get_cell_volume() == np.float(163.22171870360754)


@pytest.mark.parametrize('vasprun_parser', [('basic', {})], indirect=True)
def test_final_force(fresh_aiida_env, vasprun_parser):
    """Test that the forces are returned correctly."""

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['forces'])
    data_obj = NodeComposer.compose('array', inputs)
    # check object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    forces_check = np.array([[-0.24286901, 0., 0.], [-0.24286901, 0., 0.], [3.41460162, 0., 0.], [0.44305748, 0., 0.],
                             [-0.73887169, 0.43727184, 0.43727184], [-0.94708885, -0.85011586, 0.85011586],
                             [-0.94708885, 0.85011586, -0.85011586], [-0.73887169, -0.43727184, -0.43727184]])

    forces = data_obj.get_array('final')
    # check first, third and last position
    np.testing.assert_allclose(forces[0], forces_check[0], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(forces[2], forces_check[2], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(forces[7], forces_check[7], atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('basic', {})], indirect=True)
def test_final_stress(fresh_aiida_env, vasprun_parser):
    """Test that the stress are returned correctly."""

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['stress'])
    data_obj = NodeComposer.compose('array', inputs)
    # check object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    stress_check = np.array([[-0.38703740, 0.00000000, 0.00000000], [0.00000000, 12.52362644, -25.93894358],
                             [0.00000000, -25.93894358, 12.52362644]])
    stress = data_obj.get_array('final')
    # check entries
    np.testing.assert_allclose(stress[0], stress_check[0], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(stress[1], stress_check[1], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(stress[2], stress_check[2], atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('basic', {})], indirect=True)
def test_traj_forces(fresh_aiida_env, vasprun_parser):
    """
    Check that the parsed forces in TrajectoryData are of type ArrayData.

    Also check that the entries are as expected, e.g. correct value and
    that the first and last entry is the same (static run).

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['trajectory'])
    data_obj = NodeComposer.compose('array.trajectory', inputs)

    # test object
    ref_obj = get_data_class('array.trajectory')
    assert isinstance(data_obj, ref_obj)
    data_obj_arr = data_obj.get_array('forces')
    # test entries
    np.testing.assert_allclose(data_obj_arr[0][0], np.array([-0.24286901, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj_arr[0][-1], np.array([-0.73887169, -0.43727184, -0.43727184]), atol=0., rtol=1.0e-7)
    # test object
    ref_obj = get_data_class('array.trajectory')
    assert isinstance(data_obj, ref_obj)
    data_obj = data_obj.get_array('forces')
    # test entries
    np.testing.assert_allclose(data_obj[0][0], np.array([-0.24286901, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj[0][-1], np.array([-0.73887169, -0.43727184, -0.43727184]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj[0][-1], data_obj[1][-1], atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj[0][0], data_obj[1][0], atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('relax', {})], indirect=True)
def test_traj_forces_result_relax(fresh_aiida_env, vasprun_parser):
    """
    Check that the parsed forces in TrajectoryData are of type ArrayData.

    Also check that the entries are as expected, e.g. correct value and
    that the first and last entry is the same (static run).

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['trajectory'])
    data_obj = NodeComposer.compose('array.trajectory', inputs)
    # test object
    ref_obj = get_data_class('array.trajectory')
    assert isinstance(data_obj, ref_obj)
    data_obj_arr = data_obj.get_array('forces')
    # test shape of array
    assert data_obj_arr.shape == (19, 8, 3)
    # test a few entries (first and last atom)
    np.testing.assert_allclose(data_obj_arr[0][0], np.array([-2.42632080e-01, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj_arr[0][-1], np.array([-7.38879520e-01, -4.37063010e-01, -4.37063010e-01]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj_arr[-1][0], np.array([1.55852000e-03, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj_arr[-1][-1], np.array([-1.75970000e-03, 1.12150000e-04, 1.12150000e-04]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('relax', {})], indirect=True)
def test_unitcells_result_relax(fresh_aiida_env, vasprun_parser):
    """
    Check that the parsed unitcells are of type ArrayData.

    Also check that the entries are as expected, e.g. correct value and
    that the first and last entry is the same (static run).

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['trajectory'])
    data_obj = NodeComposer.compose('array.trajectory', inputs)
    # test object
    ref_obj = get_data_class('array.trajectory')
    assert isinstance(data_obj, ref_obj)
    data_obj_arr = data_obj.get_array('cells')
    # test shape of array
    assert data_obj_arr.shape == (19, 3, 3)
    # test a few entries (first and last vector)
    np.testing.assert_allclose(data_obj_arr[0][0], np.array([5.46503124e+00, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj_arr[0][-1], np.array([0.0, 0.0, 5.46503124e+00]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj_arr[-1][0], np.array([5.46702248e+00, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj_arr[-1][-1], np.array([0.0, 2.19104000e-03, 5.46705225e+00]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('relax', {})], indirect=True)
def test_positions_result_relax(fresh_aiida_env, vasprun_parser):
    """
    Check that the parsed positions are of type ArrayData.

    Also check that the entries are as expected, e.g. correct value and
    that the first and last entry is the same (static run).

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['trajectory'])
    data_obj = NodeComposer.compose('array.trajectory', inputs)
    # test object
    ref_obj = get_data_class('array.trajectory')
    assert isinstance(data_obj, ref_obj)
    data_obj_arr = data_obj.get_array('positions')
    # test shape of array
    assert data_obj_arr.shape == (19, 8, 3)
    # test a few entries (first and last atom)
    np.testing.assert_allclose(data_obj_arr[0][0], np.array([0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj_arr[0][-1], np.array([0.75, 0.75, 0.25]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj_arr[-1][0], np.array([-0.00621692, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(data_obj_arr[-1][-1], np.array([0.7437189, 0.74989833, 0.24989833]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('dielectric', {})], indirect=True)
def test_dielectrics(fresh_aiida_env, vasprun_parser):
    """
    Check that the parsed dielectrics are of type ArrayData.

    Also check that the entries are as expected.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['dielectrics'])
    data_obj = NodeComposer.compose('array', inputs)
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    imag = data_obj.get_array('idiel')
    real = data_obj.get_array('rdiel')
    energy = data_obj.get_array('ediel')
    # test shape of arrays
    assert imag.shape == (1000, 6)
    assert real.shape == (1000, 6)
    assert energy.shape == (1000,)
    # test a few entries
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


@pytest.mark.parametrize('vasprun_parser', [('disp_details', {})], indirect=True)
def test_epsilon(fresh_aiida_env, vasprun_parser):
    """
    Check that epsilon is returned inside the dielectric node.

    Also check that the entries are as expected.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['dielectrics'])
    data_obj = NodeComposer.compose('array', inputs)
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    epsilon = data_obj.get_array('epsilon')
    epsilon_ion = data_obj.get_array('epsilon_ion')
    # test shape of arrays
    assert epsilon.shape == (3, 3)
    assert epsilon_ion.shape == (3, 3)
    # test a few entries
    test = np.array([[13.05544887, -0., 0.], [-0., 13.05544887, -0.], [0., 0., 13.05544887]])
    np.testing.assert_allclose(epsilon, test, atol=0., rtol=1.0e-7)
    test = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
    np.testing.assert_allclose(epsilon_ion, test, atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('localfield', {})], indirect=True)
def test_born(fresh_aiida_env, vasprun_parser):
    """
    Check that the Born effective charges are of type ArrayData.

    Also check that the entries are as expected.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['born_charges'])
    data_obj = NodeComposer.compose('array', inputs)
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    born = data_obj.get_array('born_charges')
    # test shape of array
    assert born.shape == (8, 3, 3)
    # test a few entries
    np.testing.assert_allclose(born[0][0], np.array([6.37225000e-03, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(born[0][-1], np.array([-4.21760000e-04, -2.19570210e-01, 3.20709600e-02]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(born[4][0], np.array([1.68565200e-01, -2.92058000e-02, -2.92058000e-02]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('basic', {})], indirect=True)
def test_dos(fresh_aiida_env, vasprun_parser):
    """
    Check that the density of states are of type ArrayData.

    Also check that the entries are as expected.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['dos'])
    data_obj = NodeComposer.compose('array', inputs)
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    dos = data_obj.get_array('tdos')
    energy = data_obj.get_array('energy')
    # test shape of array
    assert dos.shape == (301,)
    assert energy.shape == (301,)
    # test a few entries
    assert dos[150] == pytest.approx(4.1296)
    assert energy[150] == pytest.approx(2.3373)


@pytest.mark.parametrize('vasprun_parser', [('spin', {})], indirect=True)
def test_dos_spin(fresh_aiida_env, vasprun_parser):
    """
    Check that the density of states are of type ArrayData.

    Also check that the entries are as expected.

    This test is for spin separated systems.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['dos'])
    data_obj = NodeComposer.compose('array', inputs)
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    dos = data_obj.get_array('tdos')
    # test shape of array
    assert dos.shape == (
        2,
        1000,
    )
    # test a few entries
    assert dos[0, 500] == pytest.approx(0.9839)
    assert dos[1, 500] == pytest.approx(0.9844)


@pytest.mark.parametrize('vasprun_parser', [('partial', {})], indirect=True)
def test_pdos(fresh_aiida_env, vasprun_parser):
    """
    Check that the density of states are of type ArrayData.

    Also check that that the entries are as expected.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['dos'])
    data_obj = NodeComposer.compose('array', inputs)
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    dos = data_obj.get_array('pdos')
    energy = data_obj.get_array('energy')
    # test shape of array
    assert dos.shape == (8, 1000, 9)
    assert energy.shape == (1000,)
    # test a few entries
    np.testing.assert_allclose(dos[3, 500], np.array([0.0770, 0.0146, 0.0109, 0.0155, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(dos[7, 500], np.array([0.0747, 0.0121, 0.0092, 0.0116, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    assert energy[500] == pytest.approx(0.01)


@pytest.mark.parametrize('vasprun_parser', [('partial', {})], indirect=True)
def test_projectors(fresh_aiida_env, vasprun_parser):
    """
    Check that the projectors are of type ArrayData.

    Also check that that the entries are as expected.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['projectors'])
    data_obj = NodeComposer.compose('array', inputs)
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    proj = data_obj.get_array('projectors')
    # test shape of array
    assert proj.shape == (8, 64, 21, 9)
    # test a few entries
    np.testing.assert_allclose(proj[0, 0, 5], np.array([0.0, 0.012, 0.0123, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(proj[7, 0, 5], np.array([0.1909, 0.0001, 0.0001, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(proj[4, 3, 5], np.array([0.2033, 0.0001, 0.0001, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0]), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('basic', {})], indirect=True)
def test_bands(fresh_aiida_env, vasprun_parser):
    """
    Check that the eigenvalues are of type BandData.

    Also check that the entries are as expected including the
    occupancies.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['eigenvalues', 'kpoints', 'occupancies'])
    data_obj = NodeComposer.compose('array.bands', inputs)
    # test object
    ref_obj = get_data_class('array.bands')
    assert isinstance(data_obj, ref_obj)
    eigenocc = data_obj.get_bands(also_occupations=True)
    eigen = eigenocc[0]
    occ = eigenocc[1]
    # test shape of array
    assert eigen.shape == (1, 64, 21)
    assert occ.shape == (1, 64, 21)
    # test a few entries
    assert eigen[0, 0, 0] == pytest.approx(-6.2348)
    assert eigen[0, 0, 15] == pytest.approx(5.8956)
    assert eigen[0, 6, 4] == pytest.approx(-1.7424)
    assert occ[0, 0, 0] == pytest.approx(1.0)
    assert occ[0, 0, 15] == pytest.approx(0.6949)
    assert occ[0, 6, 4] == pytest.approx(1.0)


@pytest.mark.parametrize('vasprun_parser', [('spin', {})], indirect=True)
def test_band_properties_result(fresh_aiida_env, vasprun_parser):
    """Test the result of band_properties"""

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['band_properties'])
    data = NodeComposer.compose('dict', inputs).get_dict()['band_properties']
    assert data['cbm'] == pytest.approx(6.5536)
    assert data['vbm'] == pytest.approx(6.5105)
    assert data['is_direct_gap'] is False
    assert data['band_gap'] == pytest.approx(0.04310, rel=1e-3)


@pytest.mark.parametrize('vasprun_parser', [('spin', {})], indirect=True)
def test_eigenocc_spin_result(fresh_aiida_env, vasprun_parser):
    """
    Check that the eigenvalues are of type BandData.

    Also check that the entries are as expected, including the
    occupancies. This test is for spin separated systems.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['eigenvalues', 'kpoints', 'occupancies'])
    data_obj = NodeComposer.compose('array.bands', inputs)
    # test object
    ref_obj = get_data_class('array.bands')
    assert isinstance(data_obj, ref_obj)
    eigenocc = data_obj.get_bands(also_occupations=True)
    eigen = eigenocc[0]
    occ = eigenocc[1]
    # test shape of array
    assert eigen.shape == (2, 64, 25)
    assert occ.shape == (2, 64, 25)
    # test a few entries
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


@pytest.mark.parametrize('vasprun_parser', [('basic', {})], indirect=True)
def test_toten(fresh_aiida_env, vasprun_parser):
    """
    Check that the total energy are of type ArrayData.

    Also check that the entries are as expected.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['energies'])
    data_obj = NodeComposer.compose('array', inputs)
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    # Test that the default arrays are present
    assert set(data_obj.get_arraynames()) == set(['energy_extrapolated_electronic', 'energy_extrapolated', 'electronic_steps'])
    energies = data_obj.get_array('energy_extrapolated_electronic')
    test_array = np.array([-42.91113621])
    np.testing.assert_allclose(test_array, energies, atol=0., rtol=1.0e-7)
    # Test number of entries
    assert energies.shape == (1,)
    # Electronic steps should be one
    test_array = np.array([1])
    np.testing.assert_allclose(test_array, data_obj.get_array('electronic_steps'), atol=0., rtol=1.0e-7)
    # Testing on VASP 5 so final total energy should not be the same as the last electronic step total energy.
    test_array = np.array([-42.91113621])
    np.testing.assert_allclose(test_array, data_obj.get_array('energy_extrapolated'), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('basic', {'energy_type': ['energy_free', 'energy_no_entropy']})], indirect=True)
def test_toten_multiple(fresh_aiida_env, vasprun_parser):
    """
    Check that the total energy are of type ArrayData and that we can extract multiple total energies.

    Also check that the entries are as expected.

    """
    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['energies'])
    data_obj = NodeComposer.compose('array', inputs)
    # Test that the object is of the right type
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    assert set(data_obj.get_arraynames()) == set(
        ['electronic_steps', 'energy_free', 'energy_free_electronic', 'energy_no_entropy', 'energy_no_entropy_electronic'])
    test_array = np.array([-42.91231976])
    np.testing.assert_allclose(test_array, data_obj.get_array('energy_free'), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(test_array, data_obj.get_array('energy_free_electronic'), atol=0., rtol=1.0e-7)
    test_array = np.array([-42.90995265])
    np.testing.assert_allclose(test_array, data_obj.get_array('energy_no_entropy_electronic'), atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(test_array, data_obj.get_array('energy_no_entropy'), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('basic', {'electronic_step_energies': True})], indirect=True)
def test_toten_electronic(fresh_aiida_env, vasprun_parser):
    """
    Check that the total energy are of type ArrayData and that we have entries per electronic step

    Also check that the entries are as expected.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['energies'])
    data_obj = NodeComposer.compose('array', inputs)
    # Test that the object is of the right type
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    # Test that the default arrays are present
    assert set(data_obj.get_arraynames()) == set(['energy_extrapolated_electronic', 'energy_extrapolated', 'electronic_steps'])
    energies = data_obj.get_array('energy_extrapolated_electronic')
    test_array = np.array([-42.91113666, -42.91113621])
    np.testing.assert_allclose(test_array, energies, atol=0., rtol=1.0e-7)
    # Test number of entries
    assert energies.shape == (2,)
    # Electronic steps should be two
    test_array = np.array([2])
    np.testing.assert_allclose(test_array, data_obj.get_array('electronic_steps'), atol=0., rtol=1.0e-7)
    # Testing on VASP 5 so final total energy should not be the same as the last electronic step total energy.
    test_array = np.array([-42.91113621])
    np.testing.assert_allclose(test_array, data_obj.get_array('energy_extrapolated'), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('relax', {})], indirect=True)
def test_toten_relax(fresh_aiida_env, vasprun_parser):
    """
    Check that the total energies are of type ArrayData for runs with multiple ionic steps.

    Also check that the entries are as expected.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['energies'])
    data_obj = NodeComposer.compose('array', inputs)
    # Test that the object is of the right type
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    assert set(data_obj.get_arraynames()) == set(['energy_extrapolated_electronic', 'energy_extrapolated', 'electronic_steps'])
    energies = data_obj.get_array('energy_extrapolated_electronic')
    test_array = np.array([
        -42.91113348, -43.27757545, -43.36648855, -43.37734069, -43.38062479, -43.38334165, -43.38753003, -43.38708193, -43.38641449,
        -43.38701639, -43.38699488, -43.38773717, -43.38988315, -43.3898822, -43.39011239, -43.39020751, -43.39034244, -43.39044584,
        -43.39087657
    ])
    # Test energies
    np.testing.assert_allclose(test_array, energies, atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(test_array, data_obj.get_array('energy_extrapolated'), atol=0., rtol=1.0e-7)
    # Test number of entries
    assert energies.shape == test_array.shape
    # Electronic steps should be entries times one
    np.testing.assert_allclose(np.ones(19, dtype=int), data_obj.get_array('electronic_steps'), atol=0., rtol=1.0e-7)
    # Testing on VASP 5 so final total energy should not be the same as the last electronic step total energy.


@pytest.mark.parametrize('vasprun_parser', [('relax', {'electronic_step_energies': True})], indirect=True)
def test_toten_relax_electronic(fresh_aiida_env, vasprun_parser):
    """
    Check that the total energies are of type ArrayData for runs with multiple ionic steps and electronic steps.

    Also check that the entries are as expected.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['energies'])
    data_obj = NodeComposer.compose('array', inputs)
    # Test that the object is of the right type
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    assert set(data_obj.get_arraynames()) == set(['energy_extrapolated_electronic', 'energy_extrapolated', 'electronic_steps'])
    energies = data_obj.get_array('energy_extrapolated_electronic')
    test_array_energies_electronic = [
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
    for ionic_step in test_array_energies_electronic:
        test_array_energies_flattened = np.append(test_array_energies_flattened, ionic_step)
    assert energies.shape == test_array_energies_flattened.shape
    np.testing.assert_allclose(test_array_energies_flattened, energies, atol=0., rtol=1.0e-7)
    np.testing.assert_allclose(test_array_steps, data_obj.get_array('electronic_steps'), atol=0., rtol=1.0e-7)
    test_array_energies = np.array([x[-1] for x in test_array_energies_electronic])
    # Testing on VASP 5 so final total energy should not be the same as the last electronic step total energy.
    np.testing.assert_allclose(test_array_energies, data_obj.get_array('energy_extrapolated'), atol=0., rtol=1.0e-7)


@pytest.mark.parametrize('vasprun_parser', [('disp', {})], indirect=True)
def test_hessian(fresh_aiida_env, vasprun_parser):
    """
    Check that the Hessian matrix are of type ArrayData.

    Also check that the entries are as expected.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['hessian'])
    data_obj = NodeComposer.compose('array', inputs)
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    hessian = data_obj.get_array('hessian')
    # test shape
    assert hessian.shape == (24, 24)
    # test a few entries
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


@pytest.mark.parametrize('vasprun_parser', [('disp', {})], indirect=True)
def test_dynmat(fresh_aiida_env, vasprun_parser):
    """
    Check parsing of the dynamical eigenvectors and eigenvalues.

    That it is of type ArrayData and that the entries are as expected.

    """

    inputs = get_node_composer_inputs_from_file_parser(vasprun_parser, quantity_keys=['dynmat'])
    data_obj = NodeComposer.compose('array', inputs)
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
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
