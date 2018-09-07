"""
Test submitting a ConvergenceWorkChain.

Only `run` currently works.

"""
# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name, too-many-statements
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
from aiida_vasp.utils.aiida_utils import create_authinfo


@pytest.mark.wf
@pytest.mark.skipif(aiida_version() < cmp_version('1.0.0a1'), reason='work.Runner not available before 1.0.0a1')
def test_converge_wc(fresh_aiida_env, vasp_params, potentials, mock_vasp):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import WorkflowFactory, Code
    from aiida import work

    rmq_config = None
    runner = work.Runner(poll_interval=0., rmq_config=rmq_config, enable_persistence=True)
    work.set_runner(runner)

    workchain = WorkflowFactory('vasp.converge')

    mock_vasp.store()
    comp = mock_vasp.get_computer()
    create_authinfo(computer=comp).store()

    structure = PoscarParser(file_path=data_path('test_relax_wc', 'inp', 'POSCAR')).get_quantity('poscar-structure', {})['poscar-structure']
    incar = IncarParser(file_path=data_path('test_relax_wc', 'inp', 'INCAR')).get_quantity('incar', {})['incar'].get_dict()
    incar = {k: v for k, v in incar.items() if k not in ['isif', 'ibrion']}
    incar['system'] = 'test-case:test_relax_wc'

    restart_clean_workdir = get_data_node('bool', False)
    restart_clean_workdir.store()

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    if 'encut' in incar:
        del incar['encut']
    inputs.incar = get_data_node('parameter', dict=incar)
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
    assert 'output_convergence_data' in results
    assert 'output_structure_relaxed' in results
    conv_data = results['output_convergence_data']
    try:
        conv_data.get_array('pw_regular')
    except KeyError:
        pytest.fail('Did not find pw_regular in output_convergence_data')
    try:
        conv_data.get_array('kpoints_regular')
    except KeyError:
        pytest.fail('Did not find kpoints_regular in output_convergence_data')
    try:
        conv_data.get_array('pw_displacement')
    except KeyError:
        pytest.fail('Did not find pw_displacement in output_convergence_data')
    try:
        conv_data.get_array('kpoints_displacement')
    except KeyError:
        pytest.fail('Did not find kpoints_displacement in output_convergence_data')
    try:
        conv_data.get_array('pw_compression')
    except KeyError:
        pytest.fail('Did not find pw_compression in output_convergence_data')
    try:
        conv_data.get_array('kpoints_compression')
    except KeyError:
        pytest.fail('Did not find kpoints_compression in output_convergence_data')
