"""Test the vasprun.xml io interface"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy as np

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_class


def test_parse_vasprun(vasprun_parser):
    """Parse a reference vasprun.xml file with the VasprunParser and compare the result to a reference string."""

    quantity = vasprun_parser.get_quantity('occupations')
    data_obj = quantity['occupations']
    occ = data_obj.get_array('total')
    occupations = np.array([[[1., 1., 1., 1., 0.6667, 0.6667, 0.6667, -0., -0., -0.]]])
    assert occ.all() == occupations.all()
    # eFL: How do we want to store scalar values?
    #assert  == 7.29482275


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_parameter_results(vasprun_parser):
    """
    Test that the parameter node is a ParametersData instance.

    Should contain the Fermi level.

    """

    vasprun_parser.settings.update_with({'output_params': ['fermi_level', 'total_energies', 'maximum_force', 'maximum_stress']})
    quantity = vasprun_parser.get_quantity('parameters')
    data_obj = quantity['parameters']
    ref_class = get_data_class('parameter')
    assert isinstance(data_obj, ref_class)
    data_dict = data_obj.get_dict()
    assert data_dict['fermi_level'] == 5.96764939
    assert data_dict['total_energies']['energy_no_entropy'] == -42.91113621
    assert data_dict['maximum_stress'] == 28.803993008871014
    assert data_dict['maximum_force'] == 3.41460162


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_kpoints_result(vasprun_parser):
    """Test that the kpoints result node is a KpointsData instance."""

    quantity = vasprun_parser.get_quantity('kpoints')
    data_obj = quantity['kpoints']
    ref_class = get_data_class('array.kpoints')
    assert isinstance(data_obj, ref_class)
    assert np.all(data_obj.get_kpoints()[0] == np.array([0.0, 0.0, 0.0]))
    assert np.all(data_obj.get_kpoints()[-1] == np.array([0.42857143, -0.42857143, 0.42857143]))


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_structure_result(vasprun_parser):
    """
    Test that the structure result node is a StructureData instance.

    Also check various other important properties.

    """

    quantity = vasprun_parser.get_quantity('structure')
    data_obj = quantity['structure']
    # check object
    ref_obj = get_data_class('structure')
    assert isinstance(data_obj, ref_obj)
    # check cell
    unitcell = data_obj.cell
    assert np.all(unitcell[0] == np.array([5.46503124, 0.0, 0.0]))
    assert np.all(unitcell[1] == np.array([0.0, 5.46503124, 0.0]))
    assert np.all(unitcell[2] == np.array([0.0, 0.0, 5.46503124]))
    # check first and last position
    assert np.all(data_obj.sites[0].position == np.array([0.0, 0.0, 0.0]))
    assert np.all(data_obj.sites[7].position == np.array([4.09877343, 4.09877343, 1.36625781]))
    # check volume
    assert data_obj.get_cell_volume() == np.float(163.22171870360754)


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_forces_result(vasprun_parser):
    """
    Check that the parsed forces are of type ArrayData.

    Also check that the entries are as expected, e.g. correct value and
    that the first and last entry is the same (static run).

    """

    quantity = vasprun_parser.get_quantity('trajectory')
    data_obj, data_obj_arr = quantity['trajectory']
    # test object (a tupple as we store trajectory as an array as well)
    ref_obj = get_data_class('array')
    assert isinstance(data_obj_arr, ref_obj)
    data_obj_arr = data_obj_arr.get_array('forces')
    # test entries
    assert np.all(data_obj_arr[0][0] == np.array([-0.24286901, 0.0, 0.0]))
    assert np.all(data_obj_arr[0][-1] == np.array([-0.73887169, -0.43727184, -0.43727184]))
    # test object
    ref_obj = get_data_class('array.trajectory')
    assert isinstance(data_obj, ref_obj)
    data_obj = data_obj.get_array('forces')
    # test entries
    assert np.all(data_obj[0][0] == np.array([-0.24286901, 0.0, 0.0]))
    assert np.all(data_obj[0][-1] == np.array([-0.73887169, -0.43727184, -0.43727184]))
    assert np.all(data_obj[0][-1] == data_obj[1][-1])
    assert np.all(data_obj[0][0] == data_obj[1][0])


@pytest.mark.parametrize(['vasprun_parser'], [('relax',)], indirect=True)
def test_forces_result_relax(vasprun_parser):
    """
    Check that the parsed forces are of type ArrayData.

    Also check that the entries are as expected, e.g. correct value and
    that the first and last entry is the same (static run).

    """

    quantity = vasprun_parser.get_quantity('trajectory')
    _, data_obj_arr = quantity['trajectory']
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj_arr, ref_obj)
    data_obj_arr = data_obj_arr.get_array('forces')
    # test shape of array
    assert data_obj_arr.shape == (19, 8, 3)
    # test a few entries (first and last atom)
    assert np.all(data_obj_arr[0][0] == np.array([-2.42632080e-01, 0.0, 0.0]))
    assert np.all(data_obj_arr[0][-1] == np.array([-7.38879520e-01, -4.37063010e-01, -4.37063010e-01]))
    assert np.all(data_obj_arr[-1][0] == np.array([1.55852000e-03, 0.0, 0.0]))
    assert np.all(data_obj_arr[-1][-1] == np.array([-1.75970000e-03, 1.12150000e-04, 1.12150000e-04]))


@pytest.mark.parametrize(['vasprun_parser'], [('relax',)], indirect=True)
def test_unitcells_result_relax(vasprun_parser):
    """
    Check that the parsed unitcells are of type ArrayData.

    Also check that the entries are as expected, e.g. correct value and
    that the first and last entry is the same (static run).

    """

    quantity = vasprun_parser.get_quantity('trajectory')
    _, data_obj_arr = quantity['trajectory']
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj_arr, ref_obj)
    data_obj_arr = data_obj_arr.get_array('cells')
    # test shape of array
    assert data_obj_arr.shape == (19, 3, 3)
    # test a few entries (first and last vector)
    assert np.all(data_obj_arr[0][0] == np.array([5.46503124e+00, 0.0, 0.0]))
    assert np.all(data_obj_arr[0][-1] == np.array([0.0, 0.0, 5.46503124e+00]))
    assert np.all(data_obj_arr[-1][0] == np.array([5.46702248e+00, 0.0, 0.0]))
    assert np.all(data_obj_arr[-1][-1] == np.array([0.0, 2.19104000e-03, 5.46705225e+00]))


@pytest.mark.parametrize(['vasprun_parser'], [('relax',)], indirect=True)
def test_positions_result_relax(vasprun_parser):
    """
    Check that the parsed positions are of type ArrayData.

    Also check that the entries are as expected, e.g. correct value and
    that the first and last entry is the same (static run).

    """

    quantity = vasprun_parser.get_quantity('trajectory')
    _, data_obj_arr = quantity['trajectory']
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj_arr, ref_obj)
    data_obj_arr = data_obj_arr.get_array('positions')
    # test shape of array
    assert data_obj_arr.shape == (19, 8, 3)
    # test a few entries (first and last atom)
    assert np.all(data_obj_arr[0][0] == np.array([0.0, 0.0, 0.0]))
    assert np.all(data_obj_arr[0][-1] == np.array([0.75, 0.75, 0.25]))
    assert np.all(data_obj_arr[-1][0] == np.array([-0.00621692, 0.0, 0.0]))
    assert np.all(data_obj_arr[-1][-1] == np.array([0.7437189, 0.74989833, 0.24989833]))


@pytest.mark.parametrize(['vasprun_parser'], [('dielectric',)], indirect=True)
def test_dielectrics_result(vasprun_parser):
    """
    Check that the parsed dielectrics are of type ArrayData.

    Also check that the entries are as expected.

    """

    quantity = vasprun_parser.get_quantity('dielectrics')
    data_obj = quantity['dielectrics']
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    imag = data_obj.get_array('idiel')
    real = data_obj.get_array('rdiel')
    energy = data_obj.get_array('ediel')
    # test shape of array
    assert imag.shape == (1000, 6)
    assert real.shape == (1000, 6)
    assert energy.shape == (1000,)
    # test a few entries
    assert np.all(imag[0] == np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    assert np.all(imag[500] == np.array([0.0933, 0.0924, 0.0924, 0.0, 0.0082, 0.0]))

    assert np.all(imag[999] == np.array([0.0035, 0.0035, 0.0035, 0.0, 0.0, 0.0]))
    assert np.all(real[0] == np.array([12.0757, 11.4969, 11.4969, 0.0, 0.6477, 0.0]))
    assert np.all(real[500] == np.array([-0.5237, -0.5366, -0.5366, 0.0, 0.0134, 0.0]))
    assert np.all(real[999] == np.array([6.57100000e-01, 6.55100000e-01, 6.55100000e-01, 0.0, -1.00000000e-04, 0.0]))
    assert energy[500] == 10.2933


@pytest.mark.parametrize(['vasprun_parser'], [('localfield',)], indirect=True)
def test_born_result(vasprun_parser):
    """
    Check that the Born effective charges are of type ArrayData.

    Also check that the entries are as expected.

    """

    quantity = vasprun_parser.get_quantity('born_charges')
    data_obj = quantity['born_charges']
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    born = data_obj.get_array('born_charges')
    # test shape of array
    assert born.shape == (8, 3, 3)
    # test a few entries
    assert np.all(born[0][0] == np.array([6.37225000e-03, 0.0, 0.0]))
    assert np.all(born[0][-1] == np.array([-4.21760000e-04, -2.19570210e-01, 3.20709600e-02]))
    assert np.all(born[4][0] == np.array([1.68565200e-01, -2.92058000e-02, -2.92058000e-02]))


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_dos_result(vasprun_parser):
    """
    Check that the density of states are of type ArrayData.

    Also check that the entries are as expected.

    """

    quantity = vasprun_parser.get_quantity('dos')
    data_obj = quantity['dos']
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    dos = data_obj.get_array('tdos')
    energy = data_obj.get_array('energy')
    # test shape of array
    assert dos.shape == (301,)
    assert energy.shape == (301,)
    # test a few entries
    assert dos[150] == 4.1296
    assert energy[150] == 2.3373


@pytest.mark.parametrize(['vasprun_parser'], [('spin',)], indirect=True)
def test_dos_spin_result(vasprun_parser):
    """
    Check that the density of states are of type ArrayData.

    Also check that the entries are as expected.

    This test is for spin separated systems.

    """

    quantity = vasprun_parser.get_quantity('dos')
    data_obj = quantity['dos']
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
    assert dos[0, 500] == 0.9839
    assert dos[1, 500] == 0.9844


@pytest.mark.parametrize(['vasprun_parser'], [('partial',)], indirect=True)
def test_pdos_result(vasprun_parser):
    """
    Check that the density of states are of type ArrayData.

    Also check that that the entries are as expected.

    """

    quantity = vasprun_parser.get_quantity('dos')
    data_obj = quantity['dos']
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    dos = data_obj.get_array('pdos')
    energy = data_obj.get_array('energy')
    # test shape of array
    assert dos.shape == (8, 1000, 9)
    assert energy.shape == (1000,)
    # test a few entries
    assert np.all(dos[3, 500] == np.array([0.0770, 0.0146, 0.0109, 0.0155, 0.0, 0.0, 0.0, 0.0, 0.0]))
    assert np.all(dos[7, 500] == np.array([0.0747, 0.0121, 0.0092, 0.0116, 0.0, 0.0, 0.0, 0.0, 0.0]))
    assert energy[500] == 0.01


@pytest.mark.parametrize(['vasprun_parser'], [('partial',)], indirect=True)
def test_projectors_result(vasprun_parser):
    """
    Check that the projectors are of type ArrayData.

    Also check that that the entries are as expected.

    """

    quantity = vasprun_parser.get_quantity('projectors')
    data_obj = quantity['projectors']
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    proj = data_obj.get_array('projectors')
    # test shape of array
    assert proj.shape == (8, 64, 21, 9)
    # test a few entries
    assert np.all(proj[0, 0, 5] == np.array([0.0, 0.012, 0.0123, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    assert np.all(proj[7, 0, 5] == np.array([0.1909, 0.0001, 0.0001, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0]))
    assert np.all(proj[4, 3, 5] == np.array([0.2033, 0.0001, 0.0001, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0]))


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_bands_result(vasprun_parser):
    """
    Check that the eigenvalues are of type BandData.

    Also check that the entries are as expected including the
    occupancies.

    """

    quantity = vasprun_parser.get_quantity('bands')
    data_obj = quantity['bands']
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
    assert eigen[0, 0, 0] == -6.2348
    assert eigen[0, 0, 15] == 5.8956
    assert eigen[0, 6, 4] == -1.7424
    assert occ[0, 0, 0] == 1.0
    assert occ[0, 0, 15] == 0.6949
    assert occ[0, 6, 4] == 1.0


@pytest.mark.parametrize(['vasprun_parser'], [('spin',)], indirect=True)
def test_eigenocc_spin_result(vasprun_parser):
    """
    Check that the eigenvalues are of type BandData.

    Also check that the entries are as expected, including the
    occupancies. This test is for spin separated systems.

    """

    quantity = vasprun_parser.get_quantity('bands')
    data_obj = quantity['bands']
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
    assert eigen[0, 0, 0] == -6.2363
    assert eigen[0, 0, 15] == 5.8939
    assert eigen[0, 6, 4] == -1.7438
    assert eigen[1, 0, 0] == -6.2357
    assert eigen[1, 0, 15] == 5.8946
    assert eigen[1, 6, 4] == -1.7432
    assert occ[0, 0, 0] == 1.0
    assert occ[0, 0, 15] == 0.6955
    assert occ[0, 6, 4] == 1.0
    assert occ[1, 0, 0] == 1.0
    assert occ[1, 0, 15] == 0.6938
    assert occ[1, 6, 4] == 1.0


@pytest.mark.parametrize(['vasprun_parser'], [('basic',)], indirect=True)
def test_toten_result(vasprun_parser):
    """
    Check that the total energy are of type ArrayData.

    Also check that the entries are as expected.

    """

    quantity = vasprun_parser.get_quantity('energies')
    data_obj = quantity['energies']
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    energies = data_obj.get_array('energy_no_entropy')
    # test number of entries
    assert energies.shape == (1,)
    # check energy
    assert energies[0] == -42.91113621


@pytest.mark.parametrize(['vasprun_parser'], [('relax',)], indirect=True)
def test_totens_relax_result(vasprun_parser):
    """
    Check that the total energies are of type ArrayData.

    Also check that the entries are as expected.

    """

    quantity = vasprun_parser.get_quantity('energies')
    data_obj = quantity['energies']
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    energies = data_obj.get_array('energy_no_entropy')
    # test number of entries
    assert energies.shape == (19,)
    # test a few entries
    assert energies[0] == -42.91113348
    assert energies[3] == -43.37734069
    assert energies[-1] == -43.39087657


@pytest.mark.parametrize(['vasprun_parser'], [('disp',)], indirect=True)
def test_hessian_result(vasprun_parser):
    """
    Check that the Hessian matrix are of type ArrayData.

    Also check that the entries are as expected.

    """

    quantity = vasprun_parser.get_quantity('hessian')
    data_obj = quantity['hessian']
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    hessian = data_obj.get_array('hessian')
    # test shape
    assert hessian.shape == (24, 24)
    # test a few entries
    assert np.all(hessian[0] == np.array([
        -4.63550410e-01, 0.00000000e+00, 0.00000000e+00, -5.91774100e-02, 0.00000000e+00, 0.00000000e+00, 3.09711000e-02, 0.00000000e+00,
        0.00000000e+00, 3.20435400e-02, 0.00000000e+00, 0.00000000e+00, 1.15129840e-01, -8.16138200e-02, 8.17234700e-02, 1.14879520e-01,
        8.11324800e-02, 8.27409500e-02, 1.14879520e-01, -8.11324800e-02, -8.27409500e-02, 1.15129840e-01, 8.16138200e-02, -8.17234700e-02
    ]))
    assert np.all(hessian[-2] == np.array([
        8.16138200e-02, 1.15195590e-01, -8.38411100e-02, -8.17234700e-02, 1.14875090e-01, -8.53388100e-02, 3.46686900e-02, 7.00672700e-02,
        2.54288300e-02, -8.26222700e-02, 1.16185510e-01, 7.95575600e-02, -3.05970000e-04, 3.16827300e-02, 2.86379000e-03, 5.42080000e-04,
        3.27613500e-02, 1.12576000e-03, -1.34305000e-03, -5.86811100e-02, 2.83374000e-03, 4.91688400e-02, -4.22101090e-01, 5.73736900e-02
    ]))


@pytest.mark.parametrize(['vasprun_parser'], [('disp',)], indirect=True)
def test_dynmat_result(vasprun_parser):
    """
    Check parsing of the dynamical eigenvectors and eigenvalues.

    That it is of type ArrayData and that the entries are as expected.

    """

    quantity = vasprun_parser.get_quantity('dynmat')
    data_obj = quantity['dynmat']
    # test object
    ref_obj = get_data_class('array')
    assert isinstance(data_obj, ref_obj)
    dynvec = data_obj.get_array('dynvec')
    dyneig = data_obj.get_array('dyneig')
    # test shape
    assert dynvec.shape == (24, 24)
    assert dyneig.shape == (24,)
    # test a few entries
    assert np.all(dynvec[0] == np.array([
        7.28517310e-17, 7.25431601e-02, -4.51957676e-02, 1.15412776e-16, 4.51957676e-02, -7.25431601e-02, -1.37347223e-16, 5.16257351e-01,
        -5.16257351e-01, 8.16789156e-17, 8.95098005e-02, -8.95098005e-02, -4.43838008e-17, -6.38031134e-02, 6.38031134e-02, -1.80132830e-01,
        -2.97969516e-01, 2.97969516e-01, 1.80132830e-01, -2.97969516e-01, 2.97969516e-01, -2.09989969e-16, -6.38031134e-02, 6.38031134e-02
    ]))
    assert np.all(dynvec[4] == np.array([
        -5.29825122e-13, -2.41759046e-01, -3.28913434e-01, -5.30734671e-13, -3.28913434e-01, -2.41759046e-01, 3.26325910e-13,
        -3.80807441e-02, -3.80807441e-02, -9.22956103e-13, -2.99868012e-01, -2.99868012e-01, 1.64418993e-01, 1.81002749e-01, 1.81002749e-01,
        3.11984195e-13, 2.73349550e-01, 2.73349550e-01, 2.59853610e-13, 2.73349550e-01, 2.73349550e-01, -1.64418993e-01, 1.81002749e-01,
        1.81002749e-01
    ]))
    assert dyneig[0] == -1.36621537e+00
    assert dyneig[4] == -8.48939361e-01
