"""
Test submitting a VaspWorkChain.

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
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node, get_current_user, aiida_version, cmp_version


def create_authinfo(computer):
    """
    Allow the current user to use the given computer.

    Deal with backwards compatibility down to aiida 0.11
    """
    from aiida.orm import backend as orm_backend
    authinfo = None
    if hasattr(orm_backend, 'construct_backend'):
        backend = orm_backend.construct_backend()
        authinfo = backend.authinfos.create(computer=computer, user=get_current_user())
    else:
        from aiida.backends.settings import BACKEND
        from aiida.backends.profile import BACKEND_SQLA, BACKEND_DJANGO

        if BACKEND == BACKEND_DJANGO:
            from aiida.backends.djsite.db.models import DbAuthInfo
            authinfo = DbAuthInfo(dbcomputer=computer.dbcomputer, aiidauser=get_current_user())
        elif BACKEND == BACKEND_SQLA:
            from aiida.backends.sqlalchemy.models.authinfo import DbAuthInfo
            from aiida.backends.sqlalchemy import get_scoped_session
            _ = get_scoped_session()
            authinfo = DbAuthInfo(dbcomputer=computer.dbcomputer, aiidauser=get_current_user())
    return authinfo


@pytest.mark.wf
@pytest.mark.skipif(aiida_version() < cmp_version('1.0.0a1'), reason='work.Runner not available before 1.0.0a1')
@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_base(fresh_aiida_env, vasp_params, potentials, vasp_kpoints, vasp_structure, mock_vasp):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import WorkflowFactory, Code
    from aiida import work

    rmq_config = None
    runner = work.Runner(poll_interval=0., rmq_config=rmq_config, enable_persistence=True)
    work.set_runner(runner)

    workchain = WorkflowFactory('vasp.vasp')

    mock_vasp.store()
    print(mock_vasp.get_remote_exec_path())
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
    inputs.incar = vasp_params
    inputs.kpoints = kpoints
    inputs.potcar_family = get_data_node('str', POTCAR_FAMILY_NAME)
    inputs.potcar_mapping = get_data_node('parameter', dict=POTCAR_MAP)
    inputs.max_iterations = get_data_class('int')(2)
    inputs.options = get_data_node(
        'parameter', dict={
            'queue_name': 'None',
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 1
            }
        })
    inputs.max_iterations = get_data_node('int', 1)
    inputs.settings = get_data_node('parameter', dict={'parser_settings': {'add_structure': False, 'should_parse_CONTCAR': False}})

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
