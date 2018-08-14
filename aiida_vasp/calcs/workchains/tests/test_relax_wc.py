"""
Test submitting a RelaxWorkChain.

This does not seem to work, for `submit` the daemon will not pick up the workchain
and `run` just seems to get stuck after a while.
"""
# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name
from __future__ import print_function

import pytest
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_node, get_current_user, aiida_version, cmp_version, not_ubuntu
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.poscar import PoscarParser
from aiida_vasp.io.incar import IncarParser


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
#@pytest.mark.skipif(not_ubuntu(), reason='The workchain tests currently only work on Ubuntu systems. It uses the direct scheduler.')
def test_relax_wf(fresh_aiida_env, vasp_params, potentials, mock_vasp):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import WorkflowFactory, Code
    from aiida import work

    rmq_config = None
    runner = work.Runner(poll_interval=0., rmq_config=rmq_config, enable_persistence=True)
    work.set_runner(runner)

    workchain = WorkflowFactory('vasp.relax')

    mock_vasp.store()
    print(mock_vasp.get_remote_exec_path())
    comp = mock_vasp.get_computer()
    create_authinfo(computer=comp).store()

    structure = PoscarParser(file_path=data_path('test_relax_wc', 'inp', 'POSCAR')).get_quantity('poscar-structure', {})['poscar-structure']
    kpoints = KpParser(file_path=data_path('test_relax_wc', 'inp', 'KPOINTS')).get_quantity('kpoints-kpoints', {})['kpoints-kpoints']
    incar = IncarParser(file_path=data_path('test_relax_wc', 'inp', 'INCAR')).get_quantity('incar', {})['incar'].get_dict()
    incar = {k: v for k, v in incar.items() if k not in ['isif', 'ibrion']}
    incar['system'] = 'test-case:test_relax_wc'

    restart_clean_workdir = get_data_node('bool', False)
    restart_clean_workdir.store()

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.kpoints = kpoints
    inputs.potential_family = get_data_node('str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('parameter', dict=POTCAR_MAP)
    inputs.options = get_data_node(
        'parameter', dict={
            'queue_name': 'None',
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 1
            }
        })
    inputs.restart = AttributeDict()
    inputs.verify = AttributeDict()
    inputs.restart.max_iterations = get_data_node('int', 1)
    inputs.restart.clean_workdir = restart_clean_workdir
    inputs.verify.max_iterations = get_data_node('int', 1)
    inputs.verify.clean_workdir = restart_clean_workdir
    inputs.relax = AttributeDict()
    inputs.relax.incar = get_data_node('parameter', dict=incar)
    inputs.relax.perform = get_data_node('bool', True)
    inputs.relax.convergence = AttributeDict()
    inputs.relax.convergence.shape = AttributeDict()
    inputs.relax.convergence.on = get_data_node('bool', True)
    inputs.relax.convergence.positions = get_data_node('float', 0.1)

    results = work.run(workchain, **inputs)
    assert 'output_structure_relaxed' in results
