"""Unittests for VaspCalculation"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
import contextlib
import os
import math

import pytest
from aiida.common.folders import SandboxFolder

from aiida_vasp.parsers.file_parsers.potcar import MultiPotcarIo
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC
from aiida_vasp.utils.fixtures.data import get_data_class


@ONLY_ONE_CALC
def test_write_incar(vasp_calc_and_ref):
    """Write parameters input node to INCAR file, compare to reference string."""
    vasp_calc, reference = vasp_calc_and_ref
    with managed_temp_file() as temp_file:
        vasp_calc.write_incar(temp_file)
        with open(temp_file, 'r') as result_incar_fo:
            assert result_incar_fo.read() == reference['incar']


@ONLY_ONE_CALC
def test_write_potcar(vasp_calc_and_ref):
    """Check that POTCAR is written correctly"""
    vasp_calc, _ = vasp_calc_and_ref

    with managed_temp_file() as temp_file:
        vasp_calc.write_potcar(temp_file)
        with open(temp_file, 'r') as potcar_fo:
            result_potcar = potcar_fo.read()
        assert 'In_sv' in result_potcar
        assert 'As' in result_potcar
        assert 'In_d' in result_potcar
        assert result_potcar.count('End of Dataset') == 2

        if isinstance(vasp_calc.inputs.structure, get_data_class('structure')):
            multipotcar = MultiPotcarIo.read(temp_file)
            potcar_order = [potcar.node.full_name for potcar in multipotcar.potcars]
            assert potcar_order == ['In_sv', 'As', 'In_d', 'As']


@ONLY_ONE_CALC
def test_write_chgcar(localhost_dir, vasp_calc, vasp_inputs, vasp_chgcar):
    """Test that CHGAR file is written correctly"""
    from aiida.common.folders import Folder
    chgcar, ref_chgcar = vasp_chgcar

    inputs = vasp_inputs(parameters={'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'sigma': 0.5, 'magmom': '30 * 2*0.', 'icharg': 1})

    inputs.charge_density = chgcar

    calc = vasp_calc(inputs=inputs)
    temp_folder = Folder(str(localhost_dir.dirpath()))

    calcinfo = calc.prepare_for_submission(temp_folder)

    assert 'CHGCAR' in [item[1] for item in calcinfo.local_copy_list]


@ONLY_ONE_CALC
def test_write_wavecar(localhost_dir, vasp_calc, vasp_inputs, vasp_wavecar):
    """Test that WAVECAR file is written correctly"""
    from aiida.common.folders import Folder
    wavecar, ref_wavecar = vasp_wavecar

    inputs = vasp_inputs(parameters={'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'sigma': 0.5, 'magmom': '30 * 2*0.', 'istart': 1})

    inputs.wavefunctions = wavecar

    calc = vasp_calc(inputs=inputs)

    temp_folder = Folder(str(localhost_dir.dirpath()))

    calcinfo = calc.prepare_for_submission(temp_folder)

    assert 'WAVECAR' in [item[1] for item in calcinfo.local_copy_list]


# pylint: disable=protected-access
@ONLY_ONE_CALC
def test_prepare(vasp_calc, vasp_chgcar, vasp_wavecar, vasp_inputs, localhost_dir):
    """Check that preparing creates all necessary files"""
    from aiida.common.folders import Folder
    wavecar, ref_wavecar = vasp_wavecar
    chgcar, ref_chgcar = vasp_chgcar

    inputs_dict = {'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'sima': 0.5, 'magmom': '30 * 2*0.', 'charg': 11}

    inputs = vasp_inputs(parameters=inputs_dict)
    inputs.charge_density = chgcar
    inputs.wavefunctions = wavecar

    calc = vasp_calc(inputs=inputs)
    temp_folder = Folder(str(localhost_dir.dirpath()))

    calcinfo = calc.prepare_for_submission(temp_folder)
    input_files = temp_folder.get_content_list()

    for file_name in ['INCAR', 'KPOINTS', 'POSCAR', 'POTCAR']:
        assert file_name in input_files

    assert 'EIGENVAL' in calcinfo.retrieve_list
    assert 'DOSCAR' in calcinfo.retrieve_list
    assert ('wannier90*', '.', 0) in calcinfo.retrieve_list

    inputs_dict.update({'icharg': 2})

    inputs = vasp_inputs(parameters=inputs_dict)
    inputs.charge_density = chgcar
    inputs.wavefunctions = wavecar

    calc = vasp_calc(inputs=inputs)
    temp_folder = Folder(str(localhost_dir.dirpath()))

    calcinfo = calc.prepare_for_submission(temp_folder)

    assert 'WAVECAR' in [item[1] for item in calcinfo.local_copy_list]


@ONLY_ONE_CALC
def test_verify_success(vasp_calc_and_ref):
    """Check that correct inputs are successfully verified"""
    vasp_calc, _ = vasp_calc_and_ref
    vasp_calc.verify_inputs()


@ONLY_ONE_CALC
def test_verify_fail(vasp_calc, vasp_inputs):
    """Check that incorrect inputs are not verified"""
    inputs = vasp_inputs()
    inputs.pop('kpoints')

    with pytest.raises(ValueError):
        vasp_calc(inputs=inputs)


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
