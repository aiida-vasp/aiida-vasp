"""Unittests for VaspCalculation."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import, import-outside-toplevel
import contextlib
import os

import pytest

from aiida_vasp.parsers.content_parsers.potcar import MultiPotcarIo
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node


@ONLY_ONE_CALC
def test_write_incar(vasp_calc_and_ref):
    """Write parameters input node to INCAR object, compare to reference string."""
    vasp_calc, reference = vasp_calc_and_ref
    with managed_temp_object() as temp_object:
        vasp_calc.write_incar(temp_object)
        with open(temp_object, 'r', encoding='utf8') as result_incar_fo:
            assert result_incar_fo.readlines() == reference['incar']


@ONLY_ONE_CALC
def test_write_potcar(vasp_calc_and_ref):
    """Check that POTCAR is written correctly."""
    vasp_calc, _ = vasp_calc_and_ref

    with managed_temp_object() as temp_object:
        vasp_calc.write_potcar(temp_object)
        with open(temp_object, 'r', encoding='utf8') as potcar_fo:
            result_potcar = potcar_fo.read()
        assert 'In_sv' in result_potcar
        assert 'As' in result_potcar
        assert 'In_d' in result_potcar
        assert result_potcar.count('End of Dataset') == 2

        if isinstance(vasp_calc.inputs.structure, get_data_class('core.structure')):
            multipotcar = MultiPotcarIo.read(temp_object)
            potcar_order = [potcar.node.full_name for potcar in multipotcar.potcars]
            assert potcar_order == ['In_sv', 'As', 'In_d', 'As']


@ONLY_ONE_CALC
def test_write_chgcar(localhost_dir, vasp_calc, vasp_inputs, vasp_chgcar):
    """Test that CHGAR object is written correctly."""
    from aiida.common.folders import Folder
    chgcar, _ = vasp_chgcar

    inputs = vasp_inputs(parameters={'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'sigma': 0.5, 'magmom': '30 * 2*0.', 'icharg': 1})

    inputs.charge_density = chgcar
    calc = vasp_calc(inputs=inputs)
    temp_folder = Folder(str(localhost_dir.parent))

    calcinfo = calc.prepare_for_submission(temp_folder)
    assert 'CHGCAR' in [item[1] for item in calcinfo.local_copy_list]


@ONLY_ONE_CALC
def test_write_wavecar(localhost_dir, vasp_calc, vasp_inputs, vasp_wavecar):
    """Test that WAVECAR object is written correctly."""
    from aiida.common.folders import Folder
    wavecar, _ = vasp_wavecar
    inputs = vasp_inputs(parameters={'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'sigma': 0.5, 'magmom': '30 * 2*0.', 'istart': 1})
    inputs.wavefunctions = wavecar
    calc = vasp_calc(inputs=inputs)
    temp_folder = Folder(str(localhost_dir.parent))
    calcinfo = calc.prepare_for_submission(temp_folder)

    assert 'WAVECAR' in [item[1] for item in calcinfo.local_copy_list]


@ONLY_ONE_CALC
def test_incar_validate(vasp_calc, vasp_inputs, localhost_dir):
    """Test incar with invaid tags raises exception"""
    from aiida.common import ValidationError
    from aiida.common.folders import Folder
    inputs_dict = {
        'gga': 'PE',
        'smear': 3  # <- Invalid tag
    }
    inputs = vasp_inputs(parameters=inputs_dict)
    calc = vasp_calc(inputs=inputs)

    temp_folder = Folder(str(localhost_dir.parent))
    with pytest.raises(ValidationError):
        calc.prepare_for_submission(temp_folder)


# pylint: disable=protected-access
@ONLY_ONE_CALC
def test_prepare(vasp_calc, vasp_chgcar, vasp_wavecar, vasp_inputs, localhost_dir):
    """Check that preparing creates all necessary objects."""
    from aiida.common.folders import Folder
    from aiida_vasp.calcs.vasp import VaspCalculation
    wavecar, _ = vasp_wavecar
    chgcar, _ = vasp_chgcar

    inputs_dict = {'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'sigma': 0.5, 'magmom': '30 * 2*0.', 'icharg': 11}

    inputs = vasp_inputs(parameters=inputs_dict)
    inputs.charge_density = chgcar
    inputs.wavefunctions = wavecar

    calc = vasp_calc(inputs=inputs)
    temp_folder = Folder(str(localhost_dir.parent))
    calcinfo = calc.prepare_for_submission(temp_folder)
    input_objects = temp_folder.get_content_list()

    for name in ['INCAR', 'KPOINTS', 'POSCAR', 'POTCAR']:
        assert name in input_objects

    assert 'EIGENVAL' in calcinfo.retrieve_list
    assert 'DOSCAR' in calcinfo.retrieve_list
    assert 'wannier90*' in calcinfo.retrieve_list

    assert calcinfo.codes_info[0].stdout_name == VaspCalculation._VASP_OUTPUT
    assert calcinfo.codes_info[0].join_files is True

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
def managed_temp_object():
    """Create a temp file object for a with context, delete after use."""
    import tempfile
    _, temp_object = tempfile.mkstemp()
    try:
        yield temp_object
    finally:
        os.remove(temp_object)


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
@pytest.mark.usefixtures('fresh_aiida_env')
def test_vasp_calc(run_vasp_process, aiida_caplog):
    """Test a run of a basic VASP calculation and its details."""
    from aiida_vasp.calcs.vasp import VaspCalculation
    results, node = run_vasp_process()
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
    assert 'run_status' in misc
    assert 'run_stats' in misc

    # By default we should store all always retrieve objects
    retrieve_temporary_list_ref = []
    retrieve_list_ref = VaspCalculation._ALWAYS_RETRIEVE_LIST + ['_scheduler-stdout.txt', '_scheduler-stderr.txt']
    retrieve_temporary_list = node.get_retrieve_temporary_list()
    retrieve_list = node.get_retrieve_list()
    assert retrieve_temporary_list == retrieve_temporary_list_ref
    assert set(retrieve_list_ref) == set(retrieve_list)
    objects = node.outputs.retrieved.base.repository.list_objects()
    names = [single_object.name for single_object in objects]
    # Exclude Wannier objects as they are not in the test set
    retrieve_list_ref_no_wannier = [item for item in retrieve_list_ref if 'wannier' not in item]
    assert set(names) == set(retrieve_list_ref_no_wannier)

    # Check that we always try to parse notifications
    assert misc.get('notifications') is not None

    # Check that we do not have any notifications
    assert not misc['notifications']


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_calc_delete(run_vasp_process):
    """Test a run of a basic VASP calculation where one does not want to store the always retrieved objects after parsing."""
    retrieve_list_ref = ['_scheduler-stdout.txt', '_scheduler-stderr.txt']
    inputs = {}
    inputs['settings'] = get_data_node('core.dict', dict={'ALWAYS_STORE': False})
    _, node = run_vasp_process(inputs)
    objects = node.outputs.retrieved.base.repository.list_objects()
    names = [single_object.name for single_object in objects]
    assert set(names) == set(retrieve_list_ref)


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_calc_extra(run_vasp_process):
    """Test a run of a basic VASP calculation where one wants to keep additional objects after parsing is completed."""
    # Let us add an additional object to the retrieve_list (which do not delete the object after parse)
    # and check if it is actually there
    from aiida_vasp.calcs.vasp import VaspCalculation
    inputs = {}
    extra_object_to_keep = 'POSCAR'
    inputs['settings'] = get_data_node('core.dict', dict={'ADDITIONAL_RETRIEVE_LIST': [extra_object_to_keep]})
    _, node = run_vasp_process(inputs)
    retrieve_temporary_list_ref = []
    retrieve_list_ref = VaspCalculation._ALWAYS_RETRIEVE_LIST + ['_scheduler-stdout.txt', '_scheduler-stderr.txt', 'POSCAR']
    retrieve_temporary_list = node.get_retrieve_temporary_list()
    retrieve_list = node.get_retrieve_list()
    assert retrieve_temporary_list == retrieve_temporary_list_ref
    assert set(retrieve_list_ref) == set(retrieve_list)
    objects = node.outputs.retrieved.base.repository.list_objects()
    names = [single_object.name for single_object in objects]
    # Exclude Wannier objects as they are not in the test set
    retrieve_list_ref_no_wannier = [item for item in retrieve_list_ref if 'wannier' not in item]
    assert set(names) == set(retrieve_list_ref_no_wannier)


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_calc_delete_extra(run_vasp_process):
    """Test a run of a basic VASP calculation where one wants to retrieve additional objects but not store them after parsing."""
    # Let us add an additional object to the retrieve_list (which do not delete the object after parse)
    # and check if it is actually there
    from aiida_vasp.calcs.vasp import VaspCalculation
    retrieve_list_ref = ['_scheduler-stdout.txt', '_scheduler-stderr.txt']
    inputs = {}
    extra_object_to_keep = 'POSCAR'
    inputs['settings'] = get_data_node(
        'core.dict',
        dict={
            'ALWAYS_STORE': False,
            'ADDITIONAL_RETRIEVE_TEMPORARY_LIST': [extra_object_to_keep]
        },
    )
    _, node = run_vasp_process(inputs)
    retrieve_list_ref = ['_scheduler-stdout.txt', '_scheduler-stderr.txt']
    retrieve_temporary_list_ref = VaspCalculation._ALWAYS_RETRIEVE_LIST + ['POSCAR']
    retrieve_list = node.get_retrieve_list()
    retrieve_temporary_list = node.get_retrieve_temporary_list()
    assert set(retrieve_temporary_list) == set(retrieve_temporary_list_ref)
    assert set(retrieve_list) == set(retrieve_list_ref)
    objects = node.outputs.retrieved.base.repository.list_objects()
    names = [single_object.name for single_object in objects]
    assert set(names) == set(retrieve_list_ref)


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_calc_del_str_ext(run_vasp_process):
    """Test a run of a basic VASP calculation where one wants to retrieve additional objects and store only those."""
    # Let us add an additional object to the retrieve_list (which do not delete the object after parse)
    # and check if it is actually there
    from aiida_vasp.calcs.vasp import VaspCalculation
    retrieve_list_ref = ['_scheduler-stdout.txt', '_scheduler-stderr.txt']
    inputs = {}
    extra_object_to_keep = 'POSCAR'
    inputs['settings'] = get_data_node('core.dict', dict={'ALWAYS_STORE': False, 'ADDITIONAL_RETRIEVE_LIST': [extra_object_to_keep]})
    _, node = run_vasp_process(inputs)
    retrieve_list_ref = ['_scheduler-stdout.txt', '_scheduler-stderr.txt', 'POSCAR']
    retrieve_temporary_list_ref = VaspCalculation._ALWAYS_RETRIEVE_LIST
    retrieve_list = node.get_retrieve_list()
    retrieve_temporary_list = node.get_retrieve_temporary_list()
    assert set(retrieve_temporary_list) == set(retrieve_temporary_list_ref)
    assert set(retrieve_list_ref) == set(retrieve_list)
    retrieve_list_ref = ['_scheduler-stdout.txt', '_scheduler-stderr.txt', 'POSCAR']
    objects = node.outputs.retrieved.base.repository.list_objects()
    names = [single_object.name for single_object in objects]
    assert set(names) == set(retrieve_list_ref)


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_no_potcar_in_repo(run_vasp_process):
    """Test a VASP run to verify that there is no POTCAR object in the repository."""
    # Let us add an additional object to the retrieve_list (which do not delete the object after parse)
    # and check if it is actually there
    inputs = {}
    _, node = run_vasp_process(inputs)
    repo_objects = node.base.repository.list_object_names()
    assert 'POTCAR' not in repo_objects


@pytest.mark.parametrize('test_case,expected,has_notification', [
    ['exit_codes/converged-with-error', 703, True],
    ['exit_codes/converged', 0, False],
    ['exit_codes/unfinished', 700, False],
    ['exit_codes/elec-unconverged', 701, False],
    ['exit_codes/ionic-unconverged', 702, False],
])
@pytest.mark.parametrize([
    'vasp_structure',
    'vasp_kpoints',
], [('str', 'mesh')], indirect=True)
def test_vasp_calc_exit_codes(run_vasp_process, test_case, expected, has_notification):
    """
    Test running a VASP calculation with electronic/ionic convergence problems and
    check if the exit_codes are set accordingly.
    """
    results, node = run_vasp_process(test_case=test_case)

    # Check that the standard output is there
    assert 'retrieved' in results
    assert 'misc' in results
    assert 'remote_folder' in results

    misc = results['misc'].get_dict()
    assert node.exit_status == expected
    assert bool(misc['notifications']) is has_notification


@pytest.mark.parametrize([
    'vasp_structure',
    'vasp_kpoints',
], [('str', 'mesh')], indirect=True)
def test_vasp_calc_error_suppress(run_vasp_process):
    """
    Test running a VASP calculation with electronic/ionic convergence problems and
    check if the exit_codes are set accordingly.
    """
    results, node = run_vasp_process(
        test_case='exit_codes/converged-with-error',
        settings={
            'parser_settings': {
                'critical_notifications': {
                    'add_brmix': False
                }
            },
        },
    )

    # Check that the standard output is there
    assert 'retrieved' in results
    assert 'misc' in results
    assert 'remote_folder' in results

    misc = results['misc'].get_dict()
    assert node.exit_status == 0
    assert bool(misc['notifications'])


@pytest.mark.parametrize([
    'vasp_structure',
    'vasp_kpoints',
], [('str', 'mesh')], indirect=True)
def test_vasp_calc_error_ignore_all(run_vasp_process):
    """
    Test running a VASP calculation with electronic/ionic convergence problems and
    check if the exit_codes are set accordingly.
    """
    results, node = run_vasp_process(test_case='exit_codes/converged-with-error', settings={'parser_settings': {'ignore_all_errors': True}})

    # Check that the standard output is there
    assert 'retrieved' in results
    assert 'misc' in results
    assert 'remote_folder' in results

    misc = results['misc'].get_dict()
    assert node.exit_status == 0
    assert bool(misc['notifications'])
