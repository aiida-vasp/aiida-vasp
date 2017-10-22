# pylint: disable=unused-import,unused-argument,redefined-outer-name
"""Unittests for VaspCalculation"""
import os
import contextlib

import numpy
import pytest

from aiida_vasp.utils.fixtures import aiida_env, fresh_aiida_env, localhost, vasp_params, \
    paws, vasp_structure, vasp_kpoints, vasp_code, ref_incar, localhost_dir


@pytest.fixture()
def vasp_calc_and_ref(vasp_code, vasp_params, paws, vasp_kpoints,
                      vasp_structure, ref_incar):
    """Fixture for non varying setup of a vasp calculation"""
    from aiida_vasp.calcs.vasp import VaspCalculation
    calc = VaspCalculation()
    calc.use_code(vasp_code)
    calc.set_computer(vasp_code.get_computer())
    calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine': 1})
    calc.use_parameters(vasp_params)
    calc.use_paw(paws['In'], kind='In')
    calc.use_paw(paws['As'], kind='As')
    calc.use_structure(vasp_structure)
    kpoints, ref_kpoints = vasp_kpoints
    calc.use_kpoints(kpoints)
    return calc, {'kpoints': ref_kpoints, 'incar': ref_incar}


@pytest.mark.parametrize(
    ['vasp_structure', 'vasp_kpoints'], [('cif', 'mesh'), ('str', 'list')],
    indirect=True)
def test_store(vasp_calc_and_ref):
    vasp_calc, _ = vasp_calc_and_ref
    vasp_calc.store_all()
    assert vasp_calc.pk is not None


@pytest.mark.parametrize(
    ['vasp_structure', 'vasp_kpoints'], [('cif', 'mesh')], indirect=True)
def test_write_incar(fresh_aiida_env, vasp_calc_and_ref):
    vasp_calc, reference = vasp_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    with managed_temp_file() as temp_file:
        vasp_calc.write_incar(inp, temp_file)
        with open(temp_file, 'r') as result_incar_fo:
            assert result_incar_fo.read() == reference['incar']


@pytest.mark.parametrize(
    ['vasp_structure', 'vasp_kpoints'], [('cif', 'mesh')], indirect=True)
def test_write_poscar(fresh_aiida_env, vasp_calc_and_ref):
    from ase.io.vasp import read_vasp
    vasp_calc, _ = vasp_calc_and_ref
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


def test_write_kpoints(fresh_aiida_env, vasp_calc_and_ref):
    vasp_calc, reference = vasp_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    print inp['kpoints'].get_attrs(), reference['kpoints']
    with managed_temp_file() as temp_file:
        vasp_calc.write_kpoints(inp, temp_file)
        with open(temp_file, 'r') as result_kpoints_fo:
            assert result_kpoints_fo.read() == reference['kpoints']


@pytest.mark.parametrize(
    ['vasp_structure', 'vasp_kpoints'], [('cif', 'mesh')], indirect=True)
def test_write_potcar(fresh_aiida_env, vasp_calc_and_ref):
    """Check that POTCAR is written correctly"""
    vasp_calc, _ = vasp_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    with managed_temp_file() as temp_file:
        vasp_calc.write_potcar(inp, temp_file)
        with open(temp_file, 'r') as potcar_fo:
            result_potcar = potcar_fo.read()
        assert 'In_d' in result_potcar
        assert 'As' in result_potcar
        assert result_potcar.count('End of Dataset') == 2


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
