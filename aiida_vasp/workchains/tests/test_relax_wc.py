"""
Test submitting a RelaxWorkChain.

This does not seem to work, for `submit` the daemon will not pick up the workchain
and `run` just seems to get stuck after a while.
"""
# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name,no-member
from __future__ import print_function

import pytest
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_node, aiida_version, cmp_version, create_authinfo
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.poscar import PoscarParser
from aiida_vasp.io.incar import IncarParser


@pytest.mark.wc
@pytest.mark.skipif(aiida_version() < cmp_version('1.0.0a1'), reason='work.Runner not available before 1.0.0a1')
def test_relax_wc(fresh_aiida_env, vasp_params, potentials, mock_vasp):
    #def test_relax_wc(fresh_aiida_env, vasp_params, potentials, mock_vasp, mock_relax_wc):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import WorkflowFactory, Code
    from aiida import work

    rmq_config = None
    runner = work.Runner(poll_interval=0., rmq_config=rmq_config, enable_persistence=True)
    work.set_runner(runner)

    workchain = WorkflowFactory('vasp.relax')
    #workchain = mock_relax_wc

    mock_vasp.store()
    comp = mock_vasp.get_computer()
    create_authinfo(computer=comp).store()

    structure = PoscarParser(file_path=data_path('test_relax_wc', 'inp', 'POSCAR')).structure
    kpoints = KpParser(file_path=data_path('test_relax_wc', 'inp', 'KPOINTS')).kpoints
    parameters = IncarParser(file_path=data_path('test_relax_wc', 'inp', 'INCAR')).incar
    parameters = {k: v for k, v in parameters.items() if k not in ['isif', 'ibrion']}
    parameters['system'] = 'test-case:test_relax_wc'

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.kpoints = kpoints
    inputs.potential_family = get_data_node('str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('dict', dict=POTCAR_MAP)
    inputs.options = get_data_node(
        'dict',
        dict={
            'withmpi': False,
            'queue_name': 'None',
            'max_wallclock_seconds': 1,
            'import_sys_environment': True,
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 1
            },
        })
    inputs.max_iterations = get_data_node('int', 1)
    inputs.clean_workdir = get_data_node('bool', False)
    inputs.relax_parameters = get_data_node('dict', dict=parameters)
    inputs.relax = get_data_node('bool', True)
    inputs.convergence_on = get_data_node('bool', True)
    inputs.convergence_positions = get_data_node('float', 0.1)

    results = work.run(workchain, **inputs)
    assert 'output_structure_relaxed' in results
