"""Unittests for VaspCalculation."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import, import-outside-toplevel
import contextlib
import os
import math

import pytest
from aiida.common.folders import SandboxFolder
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.parsers.file_parsers.potcar import MultiPotcarIo
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC, calc_with_retrieved
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node, create_authinfo


@ONLY_ONE_CALC
def test_write_incar(vasp_calc_and_ref):
    """Write parameters input node to INCAR file, compare to reference string."""
    vasp_calc, reference = vasp_calc_and_ref
    with managed_temp_file() as temp_file:
        vasp_calc.write_incar(temp_file)
        with open(temp_file, 'r') as result_incar_fo:
            assert result_incar_fo.readlines() == reference['vasp']


@ONLY_ONE_CALC
def test_write_potcar(vasp_calc_and_ref):
    """Check that POTCAR is written correctly."""
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
    """Test that CHGAR file is written correctly."""
    from aiida.common.folders import Folder
    chgcar, _ = vasp_chgcar

    inputs = vasp_inputs(
        parameters={'vasp': {
            'gga': 'PE',
            'gga_compat': False,
            'lorbit': 11,
            'sigma': 0.5,
            'magmom': '30 * 2*0.',
            'icharg': 1
        }})

    inputs.charge_density = chgcar
    calc = vasp_calc(inputs=inputs)
    temp_folder = Folder(str(localhost_dir.parent))

    calcinfo = calc.prepare_for_submission(temp_folder)

    assert 'CHGCAR' in [item[1] for item in calcinfo.local_copy_list]


@ONLY_ONE_CALC
def test_write_wavecar(localhost_dir, vasp_calc, vasp_inputs, vasp_wavecar):
    """Test that WAVECAR file is written correctly."""
    from aiida.common.folders import Folder
    wavecar, _ = vasp_wavecar
    inputs = vasp_inputs(
        parameters={'vasp': {
            'gga': 'PE',
            'gga_compat': False,
            'lorbit': 11,
            'sigma': 0.5,
            'magmom': '30 * 2*0.',
            'istart': 1
        }})
    inputs.wavefunctions = wavecar
    calc = vasp_calc(inputs=inputs)
    temp_folder = Folder(str(localhost_dir.parent))
    calcinfo = calc.prepare_for_submission(temp_folder)

    assert 'WAVECAR' in [item[1] for item in calcinfo.local_copy_list]


@ONLY_ONE_CALC
def test_incar_validate(vasp_calc, vasp_inputs, localhost_dir):
    """Test incar with invaid tags raises exception"""
    from aiida.common import InputValidationError
    from aiida.common.folders import Folder
    inputs_dict = {
        'vasp': {
            'gga': 'PE',
            'smear': 3  # <- Invalid tag
        }
    }
    inputs = vasp_inputs(parameters=inputs_dict)
    calc = vasp_calc(inputs=inputs)

    temp_folder = Folder(str(localhost_dir.parent))
    with pytest.raises(InputValidationError):
        calc.prepare_for_submission(temp_folder)


# pylint: disable=protected-access
@ONLY_ONE_CALC
def test_prepare(vasp_calc, vasp_chgcar, vasp_wavecar, vasp_inputs, localhost_dir):
    """Check that preparing creates all necessary files."""
    from aiida.common.folders import Folder
    wavecar, _ = vasp_wavecar
    chgcar, _ = vasp_chgcar

    inputs_dict = {'vasp': {'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'sigma': 0.5, 'magmom': '30 * 2*0.', 'icharg': 11}}

    inputs = vasp_inputs(parameters=inputs_dict)
    inputs.charge_density = chgcar
    inputs.wavefunctions = wavecar

    calc = vasp_calc(inputs=inputs)
    temp_folder = Folder(str(localhost_dir.parent))
    calcinfo = calc.prepare_for_submission(temp_folder)
    input_files = temp_folder.get_content_list()

    for file_name in ['INCAR', 'KPOINTS', 'POSCAR', 'POTCAR']:
        assert file_name in input_files

    assert 'EIGENVAL' in calcinfo.retrieve_list
    assert 'DOSCAR' in calcinfo.retrieve_list
    assert 'wannier90*' in calcinfo.retrieve_list

    inputs_dict.update({'icharg': 2})

    inputs = vasp_inputs(parameters=inputs_dict)
    inputs.charge_density = chgcar
    inputs.wavefunctions = wavecar

    calc = vasp_calc(inputs=inputs)
    temp_folder = Folder(str(localhost_dir.parent))

    calcinfo = calc.prepare_for_submission(temp_folder)

    assert 'WAVECAR' in [item[1] for item in calcinfo.local_copy_list]


@ONLY_ONE_CALC
def test_verify_success(vasp_calc_and_ref):
    """Check that correct inputs are successfully verified."""
    vasp_calc, _ = vasp_calc_and_ref
    vasp_calc.verify_inputs()


@ONLY_ONE_CALC
def test_verify_fail(vasp_calc, vasp_inputs):
    """Check that incorrect inputs are not verified."""
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


@pytest.mark.calc
@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_calc(run_vasp_calc):
    """Test a run of a basic VASP calculation and its details."""
    from aiida_vasp.calcs.vasp import VaspCalculation
    results, node = run_vasp_calc()
    assert node.exit_status == 0

    # Check that the standard output is there
    assert 'retrieved' in results
    assert 'misc' in results
    assert 'remote_folder' in results

    # Also, we need to check that there is some minimum content in misc, as it will contain
    # None if parsing somehow failed.
    misc = results['misc'].get_dict()
    assert 'total_energies' in misc
    assert 'maximum_stress' in misc

    # By default we should store all always retrieve files
    retrieve_temporary_list_ref = []
    retrieve_list_ref = VaspCalculation._ALWAYS_RETRIEVE_LIST + ['_scheduler-stdout.txt', '_scheduler-stderr.txt']
    retrieve_temporary_list = node.get_retrieve_temporary_list()
    retrieve_list = node.get_retrieve_list()
    assert retrieve_temporary_list == retrieve_temporary_list_ref
    assert set(retrieve_list_ref) == set(retrieve_list)
    files = node.outputs.retrieved.list_objects()
    file_names = [single_file.name for single_file in files]
    # Exclude Wannier files as they are not in the test set
    retrieve_list_ref_no_wannier = [item for item in retrieve_list_ref if 'wannier' not in item]
    assert set(file_names) == set(retrieve_list_ref_no_wannier)


@pytest.mark.calc
@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_calc_delete(run_vasp_calc):
    """Test a run of a basic VASP calculation where one does not want to store the always retrieved files after parsing."""
    retrieve_list_ref = ['_scheduler-stdout.txt', '_scheduler-stderr.txt']
    inputs = {}
    inputs['settings'] = get_data_node('dict', dict={'ALWAYS_STORE': False})
    _, node = run_vasp_calc(inputs)
    files = node.outputs.retrieved.list_objects()
    file_names = [single_file.name for single_file in files]
    assert set(file_names) == set(retrieve_list_ref)


@pytest.mark.calc
@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_calc_extra(run_vasp_calc):
    """Test a run of a basic VASP calculation where one wants to keep additional files after parsing is completed."""
    # Let us add an additional file to the retrieve_list (which do not delete the file after parse)
    # and check if it is actually there
    from aiida_vasp.calcs.vasp import VaspCalculation
    inputs = {}
    extra_file_to_keep = 'POSCAR'
    inputs['settings'] = get_data_node('dict', dict={'ADDITIONAL_RETRIEVE_LIST': [extra_file_to_keep]})
    _, node = run_vasp_calc(inputs)
    retrieve_temporary_list_ref = []
    retrieve_list_ref = VaspCalculation._ALWAYS_RETRIEVE_LIST + ['_scheduler-stdout.txt', '_scheduler-stderr.txt', 'POSCAR']
    retrieve_temporary_list = node.get_retrieve_temporary_list()
    retrieve_list = node.get_retrieve_list()
    assert retrieve_temporary_list == retrieve_temporary_list_ref
    assert set(retrieve_list_ref) == set(retrieve_list)
    files = node.outputs.retrieved.list_objects()
    file_names = [single_file.name for single_file in files]
    # Exclude Wannier files as they are not in the test set
    retrieve_list_ref_no_wannier = [item for item in retrieve_list_ref if 'wannier' not in item]
    assert set(file_names) == set(retrieve_list_ref_no_wannier)


@pytest.mark.calc
@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_calc_delete_extra(run_vasp_calc):
    """Test a run of a basic VASP calculation where one wants to retrieve additional files but not store them after parsing."""
    # Let us add an additional file to the retrieve_list (which do not delete the file after parse)
    # and check if it is actually there
    from aiida_vasp.calcs.vasp import VaspCalculation
    retrieve_list_ref = ['_scheduler-stdout.txt', '_scheduler-stderr.txt']
    inputs = {}
    extra_file_to_keep = 'POSCAR'
    inputs['settings'] = get_data_node('dict', dict={'ALWAYS_STORE': False, 'ADDITIONAL_RETRIEVE_TEMPORARY_LIST': [extra_file_to_keep]})
    _, node = run_vasp_calc(inputs)
    retrieve_list_ref = ['_scheduler-stdout.txt', '_scheduler-stderr.txt']
    retrieve_temporary_list_ref = VaspCalculation._ALWAYS_RETRIEVE_LIST + ['POSCAR']
    retrieve_list = node.get_retrieve_list()
    retrieve_temporary_list = node.get_retrieve_temporary_list()
    assert set(retrieve_temporary_list) == set(retrieve_temporary_list_ref)
    assert set(retrieve_list) == set(retrieve_list_ref)
    files = node.outputs.retrieved.list_objects()
    file_names = [single_file.name for single_file in files]
    assert set(file_names) == set(retrieve_list_ref)


@pytest.mark.calc
@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_calc_del_str_ext(run_vasp_calc):
    """Test a run of a basic VASP calculation where one wants to retrieve additional files and store only those."""
    # Let us add an additional file to the retrieve_list (which do not delete the file after parse)
    # and check if it is actually there
    from aiida_vasp.calcs.vasp import VaspCalculation
    retrieve_list_ref = ['_scheduler-stdout.txt', '_scheduler-stderr.txt']
    inputs = {}
    extra_file_to_keep = 'POSCAR'
    inputs['settings'] = get_data_node('dict', dict={'ALWAYS_STORE': False, 'ADDITIONAL_RETRIEVE_LIST': [extra_file_to_keep]})
    _, node = run_vasp_calc(inputs)
    retrieve_list_ref = ['_scheduler-stdout.txt', '_scheduler-stderr.txt', 'POSCAR']
    retrieve_temporary_list_ref = VaspCalculation._ALWAYS_RETRIEVE_LIST
    retrieve_list = node.get_retrieve_list()
    retrieve_temporary_list = node.get_retrieve_temporary_list()
    assert set(retrieve_temporary_list) == set(retrieve_temporary_list_ref)
    assert set(retrieve_list_ref) == set(retrieve_list)
    retrieve_list_ref = ['_scheduler-stdout.txt', '_scheduler-stderr.txt', 'POSCAR']
    files = node.outputs.retrieved.list_objects()
    file_names = [single_file.name for single_file in files]
    assert set(file_names) == set(retrieve_list_ref)


@pytest.mark.calc
@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_no_potcar_in_repo(run_vasp_calc):
    """Test a VASP run to verify that there is no POTCAR file in the repository."""
    # Let us add an additional file to the retrieve_list (which do not delete the file after parse)
    # and check if it is actually there
    inputs = {}
    _, node = run_vasp_calc(inputs)
    repo_filenames = node.list_object_names()
    assert 'POTCAR' not in repo_filenames
