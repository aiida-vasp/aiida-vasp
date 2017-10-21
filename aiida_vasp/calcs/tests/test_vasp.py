# pylint: disable=unused-import,unused-argument,redefined-outer-name
"""Unittests for VaspCalculation"""
import os
import contextlib

import numpy
import pytest

from aiida_vasp.utils.fixtures import aiida_env, fresh_aiida_env, localhost, vasp_params, \
    paws, cif_structure, kpoints_mesh, vasp_code, aiida_structure, kpoints_list, ref_incar


@pytest.fixture()
def vasp_calc_s1(vasp_code, vasp_params, paws):
    """Fixture for non varying setup of a vasp calculation"""
    from aiida_vasp.calcs.vasp import VaspCalculation
    calc = VaspCalculation()
    calc.use_code(vasp_code)
    calc.set_computer(vasp_code.get_computer())
    calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine': 1})
    calc.use_parameters(vasp_params)
    calc.use_paw(paws['In'], kind='In')
    calc.use_paw(paws['As'], kind='As')
    return calc


@pytest.fixture()
def vasp_calc_cif(vasp_calc_s1, cif_structure, kpoints_mesh):
    """Fixture for a vasp calculation with a cif structure input and mesh kpoints"""
    vasp_calc_s1.use_structure(cif_structure)
    vasp_calc_s1.use_kpoints(kpoints_mesh)
    return vasp_calc_s1


@pytest.fixture()
def vasp_calc_str(vasp_calc_s1, aiida_structure, kpoints_list):
    """Fixture for a vasp calculation with a aiida structure input and list kpoints"""
    vasp_calc_s1.use_structure(aiida_structure)
    vasp_calc_s1.use_kpoints(kpoints_list)
    return vasp_calc_s1


@pytest.fixture(params=['cif', 'str'])
def vasp_calc(request, vasp_calc_cif, vasp_calc_str):
    if request.param == 'cif':
        return vasp_calc_cif
    return vasp_calc_str


def test_store(vasp_calc):
    vasp_calc.store_all()
    assert vasp_calc.pk is not None


def test_write_incar(vasp_calc_cif, ref_incar):
    inp = vasp_calc_cif.get_inputs_dict()
    with managed_temp_file() as temp_file:
        vasp_calc_cif.write_incar(inp, temp_file)
        with open(temp_file, 'r') as result_incar_fo:
            assert result_incar_fo.read() == ref_incar


def test_write_poscar(vasp_calc):
    from ase.io.vasp import read_vasp
    inp = vasp_calc.get_inputs_dict()
    with managed_temp_file() as temp_file:
        vasp_calc.write_poscar(inp, temp_file)
        with working_directory(temp_file):
            result_ase = read_vasp(temp_file)
            ref_ase = inp['structure'].get_ase()
            assert numpy.allclose(
                result_ase.get_cell(), ref_ase.get_cell(), atol=1e-16, rtol=0)
            assert result_ase.get_chemical_formula(
            ) == ref_ase.get_chemical_formula()


@contextlib.contextmanager
def managed_temp_file():
    import tempfile
    _, temp_file = tempfile.mkstemp()
    try:
        yield temp_file
    finally:
        os.remove(temp_file)


@contextlib.contextmanager
def working_directory(new_work_dir):
    work_dir = os.getcwd()
    try:
        if os.path.isdir(new_work_dir):
            os.chdir(new_work_dir)
        else:
            os.chdir(os.path.dirname(new_work_dir))
        yield
    finally:
        os.chdir(work_dir)
