"""Unittests for Vasp2w90Calculation"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
from __future__ import absolute_import
from __future__ import print_function
import os
import re
import tempfile

import numpy
import pytest
from py import path as py_path  # pylint: disable=no-member,no-name-in-module
from aiida.common import ValidationError
from aiida.common.folders import SandboxFolder

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC, STRUCTURE_TYPES


def normalize_contents(file_contents):
    """Remove trailing zeroes after floating point and normalize trailin newline to unix standard."""
    normalized = re.sub(r'(\d*.\d*?)0+(\s)', r'\g<1>0\2', file_contents)  # remove trailing zeroes
    if not re.match(r'\n', file_contents[-1]):  # add trailing newline if necessary
        normalized += '\n'
    return normalized


def assert_contents_equivalent(contents_a, contents_b):
    """Assert equivalence of files with floating point numbers."""
    assert normalize_contents(contents_a) == normalize_contents(contents_b)


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('cif', 'mesh'), ('str', 'list')], indirect=True)
def test_store(vasp2w90_calc_and_ref):
    vasp_calc, _ = vasp2w90_calc_and_ref
    vasp_calc.store_all()
    assert vasp_calc.pk is not None


@ONLY_ONE_CALC
def test_write_incar(fresh_aiida_env, vasp2w90_calc_and_ref):
    """Write INCAR reference file and compare to reference."""
    vasp_calc, reference = vasp2w90_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    with tempfile.NamedTemporaryFile() as temp_file:
        vasp_calc.write_incar(inp, temp_file.name)
        with open(temp_file.name, 'r') as result_incar_fo:
            assert_contents_equivalent(result_incar_fo.read(), reference['incar'])


@ONLY_ONE_CALC
def test_write_win(fresh_aiida_env, vasp2w90_calc_and_ref):
    """Write wannier90.win input file and compare to reference."""
    vasp_calc, reference = vasp2w90_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        vasp_calc.write_win(inp, temp_file.name)
        with open(temp_file.name, 'r') as result_incar_fo:
            assert result_incar_fo.read() == reference['win']


@STRUCTURE_TYPES
def test_write_poscar(fresh_aiida_env, vasp2w90_calc_and_ref, vasp_structure_poscar):
    """Write POSCAR input file and compare to reference."""
    from aiida_vasp.parsers.file_parsers.poscar import PoscarParser
    vasp_calc, _ = vasp2w90_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    with tempfile.NamedTemporaryFile() as temp_file:
        vasp_calc.write_poscar(inp, temp_file.name)
        result_structure = PoscarParser(file_path=temp_file.name).structure
        ref_structure = vasp_structure_poscar.structure
        assert result_structure.cell, ref_structure.cell
        assert result_structure.get_formula() == ref_structure.get_formula()

        ref_string = vasp_structure_poscar._parsed_object.get_string()  # pylint: disable=protected-access
        with open(temp_file.name, 'r') as poscar:
            assert_contents_equivalent(poscar.read(), ref_string)


def test_write_kpoints(fresh_aiida_env, vasp2w90_calc_and_ref):
    """Write KPOINTS file and compare to reference."""
    vasp_calc, reference = vasp2w90_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    with tempfile.NamedTemporaryFile() as temp_file:
        vasp_calc.write_kpoints(inp, temp_file.name)
        with open(temp_file.name, 'r') as result_kpoints_fo:
            assert result_kpoints_fo.read() == reference['kpoints']


@ONLY_ONE_CALC
def test_write_potcar(fresh_aiida_env, vasp2w90_calc_and_ref):
    """Check that POTCAR is written correctly"""
    vasp_calc, _ = vasp2w90_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    with tempfile.NamedTemporaryFile() as temp_file:
        vasp_calc.write_potcar(inp, temp_file.name)
        with open(temp_file.name, 'r') as potcar_fo:
            result_potcar = potcar_fo.read()
        assert 'In_d' in result_potcar
        assert 'As' in result_potcar
        assert result_potcar.count('End of Dataset') == 2


@ONLY_ONE_CALC
def test_write_chgcar(fresh_aiida_env, vasp2w90_calc_and_ref, vasp_chgcar):
    """Test that CHGCAR file is written correctly"""
    vasp_calc, _ = vasp2w90_calc_and_ref
    chgcar, ref_chgcar = vasp_chgcar
    vasp_calc.use_charge_density(chgcar)
    inp = vasp_calc.get_inputs_dict()
    with tempfile.NamedTemporaryFile() as temp_file:
        vasp_calc.write_chgcar(inp, temp_file.name)
        with open(temp_file.name, 'r') as result_chgcar_fo:
            assert result_chgcar_fo.read() == ref_chgcar


@ONLY_ONE_CALC
def test_write_wavecar(fresh_aiida_env, vasp2w90_calc_and_ref, vasp_wavecar):
    """Test that WAVECAR file is written correctly"""
    vasp_calc, _ = vasp2w90_calc_and_ref
    wavecar, ref_wavecar = vasp_wavecar
    vasp_calc.use_wavefunctions(wavecar)
    inp = vasp_calc.get_inputs_dict()
    with tempfile.NamedTemporaryFile() as temp_file:
        vasp_calc.write_wavecar(inp, temp_file.name)
        with open(temp_file.name, 'r') as result_wavecar_fo:
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
def test_parse_with_retrieved(vasp_nscf_and_ref, ref_retrieved):
    """Check that parsing is successful and creates the right output links"""
    vasp_calc, _ = vasp_nscf_and_ref
    parser = vasp_calc.get_parserclass()(vasp_calc)
    success, outputs = parser.parse_with_retrieved({'retrieved': ref_retrieved})
    outputs = dict(outputs)
    assert success
    assert 'output_bands' in outputs
    assert 'output_dos' in outputs
    assert 'output_parameters' in outputs


def test_verify_success(vasp2w90_calc_and_ref):
    """Check that correct inputs are successfully verified"""
    vasp_calc, _ = vasp2w90_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    vasp_calc.verify_inputs(inp)


def test_verify_fail(vasp2w90_calc_and_ref):
    """Check that incorrect inputs are not verified"""
    vasp_calc, _ = vasp2w90_calc_and_ref
    inp = vasp_calc.get_inputs_dict()
    inp.pop('kpoints')
    with pytest.raises(ValidationError):
        vasp_calc.verify_inputs(inp)
