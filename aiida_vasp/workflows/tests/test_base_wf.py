"""
Test submitting a VaspBaseWf workflow.

This does not seem to work, for `submit` the daemon will not pick up the workchain
and `run` just seems to get stuck after a while.
"""
# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name
from __future__ import print_function
import time
import os
import uuid
import subprocess as sp

import pytest
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.aiida_utils import get_data_node, get_current_user


@pytest.mark.xfail()
@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_base(vasp_params, potentials, vasp_kpoints, vasp_structure, mock_vasp):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import WorkflowFactory, Code
    from aiida.orm.backend import construct_backend
    from aiida import work
    from aiida.work.db_types import Str, Int

    rmq_config = None
    runner = work.Runner(poll_interval=0., rmq_config=rmq_config, enable_persistence=False)
    work.set_runner(runner)

    base_wf_proc = WorkflowFactory('vasp.base')

    mock_vasp.store()
    print(mock_vasp.get_remote_exec_path())
    comp = mock_vasp.get_computer()
    backend = construct_backend()
    backend.authinfos.create(computer=comp, user=get_current_user()).store()

    # ~ os_env = os.environ.copy()
    # ~ sp.call(['verdi', 'daemon', 'start'], env=os_env)
    # ~ print sp.check_output(['verdi', 'daemon', 'status'], env=os_env)
    # ~ print sp.check_output(['which', 'verdi'], env=os_env)

    kpoints, _ = vasp_kpoints
    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = vasp_structure
    inputs.incar = vasp_params
    inputs.kpoints = kpoints
    inputs.potcar_family = Str(POTCAR_FAMILY_NAME)
    inputs.potcar_mapping = get_data_node('parameter', dict=POTCAR_MAP)
    inputs.options = get_data_node(
        'parameter', dict={
            'queue_name': 'None',
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 1
            }
        })
    inputs.max_iterations = Int(1)
    inputs.settings = get_data_node('parameter', dict={'parser_settings': {'add_structure': False, 'should_parse_CONTCAR': False}})

    # ~ workchain = run(base_wf_proc, **inputs)
    results = work.run(base_wf_proc, **inputs)
    # ~ workchain = load_node(running.pk)
    # ~ timeout = 5
    # ~ waiting_for = 0
    # ~ while not workchain.is_terminated and waiting_for < timeout:
    # ~ time.sleep(1)
    # ~ waiting_for += 1
    assert 'retrieved' in results
    assert 'output_parameters' in results
    assert 'remote_folder' in results
    # ~ assert workchain.is_finished_ok
