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
from aiida_vasp.utils.aiida_utils import get_data_node, aiida_version, cmp_version, not_ubuntu, create_authinfo
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.poscar import PoscarParser
from aiida_vasp.io.incar import IncarParser


@pytest.mark.wc
@pytest.mark.skipif(aiida_version() < cmp_version('1.0.0a1'), reason='work.Runner not available before 1.0.0a1')
#@pytest.mark.skipif(not_ubuntu(), reason='The workchain tests currently only work on Ubuntu systems. It uses the direct scheduler.')
def test_relax_wc(fresh_aiida_env, vasp_params, potentials, mock_vasp):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import WorkflowFactory, Code
    from aiida import work

    rmq_config = None
    runner = work.Runner(poll_interval=0., rmq_config=rmq_config, enable_persistence=True)
    work.set_runner(runner)

    workchain = WorkflowFactory('vasp.relax')

    mock_vasp.store()
    comp = mock_vasp.get_computer()
    create_authinfo(computer=comp).store()

    structure = PoscarParser(file_path=data_path('test_relax_wc', 'inp', 'POSCAR')).get_quantity('poscar-structure', {})['poscar-structure']
    kpoints = KpParser(file_path=data_path('test_relax_wc', 'inp', 'KPOINTS')).get_quantity('kpoints-kpoints', {})['kpoints-kpoints']
    parameters = IncarParser(file_path=data_path('test_relax_wc', 'inp', 'INCAR')).get_quantity('incar', {})['incar'].get_dict()
    parameters = {k: v for k, v in parameters.items() if k not in ['isif', 'ibrion']}
    parameters['system'] = 'test-case:test_relax_wc'

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.kpoints = kpoints
    inputs.potential_family = get_data_node('str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('parameter', dict=POTCAR_MAP)
    inputs.options = get_data_node(
        'parameter',
        dict={
            'queue_name': 'None',
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 1
            },
            'max_wallclock_seconds': 3600
        })
    inputs.max_iterations = get_data_node('int', 1)
    inputs.clean_workdir = get_data_node('bool', False)
    inputs.relax_parameters = get_data_node('parameter', dict=parameters)
    inputs.perform = get_data_node('bool', True)
    inputs.convergence_on = get_data_node('bool', True)
    inputs.convergence_positions = get_data_node('float', 0.1)

    results = work.run(workchain, **inputs)
    assert 'output_structure_relaxed' in results
    assert False
