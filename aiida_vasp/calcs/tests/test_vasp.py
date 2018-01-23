"""Unittests for VaspCalculation"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
import contextlib
import os

import pytest
from aiida.common.exceptions import ValidationError
from aiida.common.folders import SandboxFolder

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC, STRUCTURE_TYPES


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('cif', 'mesh'), ('str', 'list')], indirect=True)
def test_store(vasp_calc_and_ref):
    vasp_calc, _ = vasp_calc_and_ref
    vasp_calc.store_all()
    assert vasp_calc.pk is not None


@ONLY_ONE_CALC
def test_write_incar(fresh_aiida_env, vasp_calc_and_ref):
    """Write parameters input node to INCAR file, compare to reference string."""
    vasp_calc, reference = vasp_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    with managed_temp_file() as temp_file:
        vasp_calc.write_incar(inp, temp_file)
        with open(temp_file, 'r') as result_incar_fo:
            assert result_incar_fo.read() == reference['incar']


@STRUCTURE_TYPES
def test_write_poscar(fresh_aiida_env, vasp_calc_and_ref, vasp_structure_poscar):
    """Write structure input node to file, compare contents to reference string."""
    from pymatgen.io.vasp.inputs import Poscar
    vasp_calc, _ = vasp_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    with managed_temp_file() as temp_file:
        vasp_calc.write_poscar(inp, temp_file)
        with working_directory(temp_file):
            result_pmg = Poscar.from_file(temp_file).structure
            ref_pmg = vasp_structure_poscar.structure
            assert result_pmg.lattice, ref_pmg.lattice
            assert result_pmg.formula == ref_pmg.formula

        with open(temp_file, 'r') as poscar:
            assert poscar.read() == vasp_structure_poscar.get_string()


def test_write_kpoints(fresh_aiida_env, vasp_calc_and_ref):
    """Write the kpoints input node to a file, compare content to a reference string."""
    vasp_calc, reference = vasp_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    print inp['kpoints'].get_attrs(), reference['kpoints']
    with managed_temp_file() as temp_file:
        vasp_calc.write_kpoints(inp, temp_file)
        with open(temp_file, 'r') as result_kpoints_fo:
            assert result_kpoints_fo.read() == reference['kpoints']


@ONLY_ONE_CALC
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


@ONLY_ONE_CALC
def test_write_chgcar(fresh_aiida_env, vasp_calc_and_ref, vasp_chgcar):
    """Test that CHGAR file is written correctly"""
    vasp_calc, _ = vasp_calc_and_ref
    chgcar, ref_chgcar = vasp_chgcar
    vasp_calc.use_charge_density(chgcar)
    inp = vasp_calc.get_inputs_dict()
    with managed_temp_file() as temp_file:
        vasp_calc.write_chgcar(inp, temp_file)
        with open(temp_file, 'r') as result_chgcar_fo:
            assert result_chgcar_fo.read() == ref_chgcar


@ONLY_ONE_CALC
def test_write_wavecar(fresh_aiida_env, vasp_calc_and_ref, vasp_wavecar):
    """Test that CHGAR file is written correctly"""
    vasp_calc, _ = vasp_calc_and_ref
    wavecar, ref_wavecar = vasp_wavecar
    vasp_calc.use_wavefunctions(wavecar)
    inp = vasp_calc.get_inputs_dict()
    with managed_temp_file() as temp_file:
        vasp_calc.write_wavecar(inp, temp_file)
        with open(temp_file, 'r') as result_wavecar_fo:
            assert result_wavecar_fo.read() == ref_wavecar


# pylint: disable=protected-access
def test_prepare(vasp_nscf_and_ref):
    """Check that preparing creates all necessary files"""
    vasp_calc, _ = vasp_nscf_and_ref
    inp = vasp_calc.get_inputs_dict()
    with SandboxFolder() as sandbox_f:
        calc_info = vasp_calc._prepare_for_submission(sandbox_f, inp)
        inputs = sandbox_f.get_content_list()
    assert set(inputs) == {'INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', 'CHGCAR', 'WAVECAR'}
    assert 'EIGENVAL' in calc_info.retrieve_list
    assert 'DOSCAR' in calc_info.retrieve_list
    assert ('wannier90*', '.', 0) in calc_info.retrieve_list

    vasp_calc.inp.parameters.update_dict({'icharg': 2})
    inp = vasp_calc.get_inputs_dict()
    with SandboxFolder() as sandbox_f:
        calc_info = vasp_calc._prepare_for_submission(sandbox_f, inp)
        inputs = sandbox_f.get_content_list()
    assert set(inputs) == {'INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', 'WAVECAR'}


@ONLY_ONE_CALC
def test_parse_with_retrieved(vasp_nscf_and_ref, ref_retrieved_nscf):
    """Check that parsing is successful and creates the right output links"""
    vasp_calc, _ = vasp_nscf_and_ref
    parser = vasp_calc.get_parserclass()(vasp_calc)
    success, outputs = parser.parse_with_retrieved({'retrieved': ref_retrieved_nscf})
    outputs = dict(outputs)
    assert success
    assert 'bands' in outputs
    assert 'dos' in outputs
    assert 'results' in outputs


def test_verify_success(vasp_calc_and_ref):
    """Check that correct inputs are successfully verified"""
    vasp_calc, _ = vasp_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    vasp_calc.verify_inputs(inp)


def test_verify_fail(vasp_calc_and_ref):
    """Check that incorrect inputs are not verified"""
    vasp_calc, _ = vasp_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    inp.pop('kpoints')
    with pytest.raises(ValidationError):
        vasp_calc.verify_inputs(inp)


@contextlib.contextmanager
def managed_temp_file():
    """Create a temp file for a with context, delete after use."""
    import tempfile
    _, temp_file = tempfile.mkstemp()
    try:
        yield temp_file
    finally:
        os.remove(temp_file)


@contextlib.contextmanager
def working_directory(new_work_dir):
    """Change working dir in a with context, change back at the end."""
    work_dir = os.getcwd()
    try:
        if os.path.isdir(new_work_dir):
            os.chdir(new_work_dir)
        else:
            os.chdir(os.path.dirname(new_work_dir))
        yield
    finally:
        os.chdir(work_dir)
