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
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_node, get_current_user, aiida_version, cmp_version
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.poscar import PoscarParser
from aiida_vasp.io.incar import IncarParser
from aiida_vasp.utils.aiida_utils import create_authinfo


@pytest.mark.wf
@pytest.mark.skipif(aiida_version() < cmp_version('1.0.0a1'), reason='work.Runner not available before 1.0.0a1')
def test_relax_wf(fresh_aiida_env, vasp_params, potentials, mock_vasp):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import WorkflowFactory, Code
    from aiida import work

    rmq_config = None
    runner = work.Runner(poll_interval=0., rmq_config=rmq_config, enable_persistence=True)
    work.set_runner(runner)

    base_wf_proc = WorkflowFactory('vasp.relax')

    mock_vasp.store()
    print(mock_vasp.get_remote_exec_path())
    comp = mock_vasp.get_computer()
    create_authinfo(computer=comp).store()

    structure = PoscarParser(file_path=data_path('test_relax_wf', 'inp', 'POSCAR')).get_quantity('poscar-structure', {})['poscar-structure']
    kpoints = KpParser(file_path=data_path('test_relax_wf', 'inp', 'KPOINTS')).get_quantity('kpoints-kpoints', {})['kpoints-kpoints']
    incar_add = IncarParser(file_path=data_path('test_relax_wf', 'inp', 'INCAR')).get_quantity('incar', {})['incar'].get_dict()
    incar_add = {k: v for k, v in incar_add.items() if k not in ['isif', 'ibrion']}
    incar_add['system'] = 'test-case:test_relax_wf'

    restart_clean_workdir = get_data_node('bool', False)
    restart_clean_workdir.store()

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.incar_add = get_data_node('parameter', dict=incar_add)
    inputs.kpoints = AttributeDict()
    inputs.kpoints.mesh = kpoints
    inputs.potcar_family = get_data_node('str', POTCAR_FAMILY_NAME)
    inputs.potcar_mapping = get_data_node('parameter', dict=POTCAR_MAP)
    inputs.options = get_data_node(
        'parameter', dict={
            'queue_name': 'None',
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 1
            }
        })
    inputs.max_iterations = get_data_node('int', 1)
    inputs.convergence = AttributeDict()
    inputs.convergence.shape = AttributeDict()
    inputs.convergence.on = get_data_node('bool', True)
    inputs.convergence.positions = get_data_node('float', 0.1)
    inputs.restart = AttributeDict()
    inputs.restart.clean_workdir = restart_clean_workdir
    inputs.relax = AttributeDict()

    results = work.run(base_wf_proc, **inputs)
    assert 'relaxed_structure' in results
