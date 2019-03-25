"""
Test submitting a VaspWorkChain.

This does not seem to work, for `submit` the daemon will not pick up the workchain
and `run` just seems to get stuck after a while.
"""
# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name
from __future__ import print_function

import pytest
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.aiida_utils import get_data_node, aiida_version, cmp_version, create_authinfo


@pytest.mark.wc
@pytest.mark.skipif(aiida_version() < cmp_version('1.0.0a1'), reason='work.Runner not available before 1.0.0a1')
@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_wc(fresh_aiida_env, vasp_params, potentials, vasp_kpoints, vasp_structure, mock_vasp):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import WorkflowFactory, Code
    from aiida import work

    rmq_config = None
    runner = work.Runner(poll_interval=0., rmq_config=rmq_config, enable_persistence=True)
    work.set_runner(runner)

    workchain = WorkflowFactory('vasp.vasp')

    mock_vasp.store()
    comp = mock_vasp.get_computer()
    create_authinfo(computer=comp, store=True)

    # ~ os_env = os.environ.copy()
    # ~ sp.call(['verdi', 'daemon', 'start'], env=os_env)
    # ~ print sp.check_output(['verdi', 'daemon', 'status'], env=os_env)
    # ~ print sp.check_output(['which', 'verdi'], env=os_env)

    kpoints, _ = vasp_kpoints
    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = vasp_structure
    inputs.parameters = vasp_params
    inputs.kpoints = kpoints
    inputs.potential_family = get_data_node('str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('dict', dict=POTCAR_MAP)
    inputs.options = get_data_node(
        'dict',
        dict={
            'withmpi': False,
            'queue_name': 'None',
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 1
            },
            'max_wallclock_seconds': 3600
        })
    inputs.settings = get_data_node('dict', dict={'parser_settings': {'add_structure': False, 'should_parse_CONTCAR': False}})
    inputs.max_iterations = get_data_node('int', 1)
    inputs.clean_workdir = get_data_node('bool', False)

    # ~ running = run(workchain, **inputs)
    running = work.run(workchain, **inputs)
    # ~ running = load_node(running.pk)
    # ~ timeout = 5
    # ~ waiting_for = 0
    # ~ while not running.is_terminated and waiting_for < timeout:
    # ~ time.sleep(1)
    # ~ waiting_for += 1
    assert 'retrieved' in running
    assert 'output_parameters' in running
    assert 'remote_folder' in running
    # ~ assert running.is_finished_ok


@pytest.mark.wc
@pytest.mark.skipif(aiida_version() < cmp_version('1.0.0a1'), reason='work.Runner not available before 1.0.0a1')
@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_wc_chgcar(fresh_aiida_env, vasp_params, potentials, vasp_kpoints, vasp_structure, mock_vasp):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import WorkflowFactory, Code
    from aiida import work

    rmq_config = None
    runner = work.Runner(poll_interval=0., rmq_config=rmq_config, enable_persistence=True)
    work.set_runner(runner)

    workchain = WorkflowFactory('vasp.vasp')

    mock_vasp.store()
    comp = mock_vasp.get_computer()
    create_authinfo(computer=comp).store()

    # ~ os_env = os.environ.copy()
    # ~ sp.call(['verdi', 'daemon', 'start'], env=os_env)
    # ~ print sp.check_output(['verdi', 'daemon', 'status'], env=os_env)
    # ~ print sp.check_output(['which', 'verdi'], env=os_env)

    kpoints, _ = vasp_kpoints
    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = vasp_structure
    inputs.parameters = vasp_params
    inputs.kpoints = kpoints
    inputs.potential_family = get_data_node('str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('dict', dict=POTCAR_MAP)
    inputs.options = get_data_node(
        'dict',
        dict={
            'withmpi': False,
            'queue_name': 'None',
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 1
            },
            'max_wallclock_seconds': 3600
        })
    inputs.settings = get_data_node('dict', dict={'ADDITIONAL_RETRIEVE_LIST': ['CHGCAR'], 'parser_settings': {'add_chgcar': True}})
    inputs.max_iterations = get_data_node('int', 1)
    inputs.clean_workdir = get_data_node('bool', False)

    running = work.run(workchain, **inputs)
    assert 'output_chgcar' in running
