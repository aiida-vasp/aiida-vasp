"""
Test submitting a ConvergenceWorkChain.

Only `run` currently works.

"""
# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name, too-many-statements, import-outside-toplevel
from __future__ import print_function

import numpy as np
import pytest
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_node
from aiida_vasp.parsers.file_parsers.kpoints import KpointsParser
from aiida_vasp.parsers.file_parsers.poscar import PoscarParser
from aiida_vasp.parsers.file_parsers.incar import IncarParser
from aiida_vasp.utils.aiida_utils import create_authinfo


def test_converge_wc(fresh_aiida_env, potentials, mock_vasp):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import Code
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run

    workchain = WorkflowFactory('vasp.converge')

    mock_vasp.store()
    create_authinfo(computer=mock_vasp.computer, store=True)

    structure = PoscarParser(file_path=data_path('test_converge_wc', 'inp', 'POSCAR')).structure
    parameters = IncarParser(file_path=data_path('test_converge_wc', 'inp', 'INCAR')).incar
    parameters['system'] = 'test-case:test_converge_wc'
    parameters = {k: v for k, v in parameters.items() if k not in ['isif', 'ibrion', 'encut', 'nsw']}

    restart_clean_workdir = get_data_node('bool', False)
    restart_clean_workdir.store()

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.parameters = get_data_node('dict', dict={'incar': parameters})
    inputs.potential_family = get_data_node('str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('dict', dict=POTCAR_MAP)
    inputs.options = get_data_node('dict',
                                   dict={
                                       'withmpi': False,
                                       'queue_name': 'None',
                                       'resources': {
                                           'num_machines': 1,
                                           'num_mpiprocs_per_machine': 1
                                       },
                                       'max_wallclock_seconds': 3600
                                   })
    inputs.max_iterations = get_data_node('int', 1)
    inputs.clean_workdir = get_data_node('bool', False)
    relax = AttributeDict()
    converge = AttributeDict()
    converge.relax = get_data_node('bool', False)
    converge.compress = get_data_node('bool', False)
    converge.displace = get_data_node('bool', False)
    converge.pwcutoff_samples = get_data_node('int', 3)
    converge.k_samples = get_data_node('int', 3)
    relax.perform = get_data_node('bool', True)
    inputs.relax = relax
    inputs.converge = converge
    inputs.verbose = get_data_node('bool', True)
    results, node = run.get_node(workchain, **inputs)
    assert node.exit_status == 0
    converge = results['converge']
    assert 'data' in converge

    conv_data = converge['data']
    try:
        'pw_regular' not in conv_data
    except KeyError:
        # pytest.fail('Did not find pw_regular in converge.data')
        pytest.xfail('Cannot find pw_regular since no test datum are supplied.')
    try:
        'kpoints_regular' not in conv_data
    except KeyError:
        # pytest.fail('Did not find kpoints_regular in converge.data')
        pytest.xfail('Cannot find kpoints_regular since no test datum are supplied.')

    assert 'pwcutoff_recommended' in converge
    try:
        _encut = converge['pwcutoff_recommended'].value
    except AttributeError:
        # pytest.fail('pwcutoff_recommended does not have the expected format')
        pytest.xfail('Cannot find pwcuoff_recommended since no test datum are supplied.')
    assert 'kpoints_recommended' in converge
    try:
        _kpoints = converge['kpoints_recommended'].get_kpoints_mesh()
    except AttributeError:
        # pytest.fail('kpoints_recommended does not have the expected format')
        pytest.xfail('Cannot find kpoints_recommended since no test datum are supplied.')


def compare_converge_datum(conv_data_actual, conv_data_expect):
    """Compare pw_data or k_data in outputs of ConvergeWorkChain with a test case."""
    assert len(conv_data_actual) == len(conv_data_expect)
    for data_actual, data_expect in zip(conv_data_actual, conv_data_expect):
        assert len(data_actual) == len(data_expect)
        for actual, expect in zip(data_actual, data_expect):
            if actual is None:
                assert expect is None
            else:
                assert np.isclose(actual, expect)


def test_converge_wc_pw(fresh_aiida_env, vasp_params, potentials, mock_vasp):
    """Test convergence workflow using mock code."""
    from aiida.orm import Code
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run

    workchain = WorkflowFactory('vasp.converge')

    mock_vasp.store()
    create_authinfo(computer=mock_vasp.computer).store()

    structure = PoscarParser(file_path=data_path('test_converge_wc/pw/200', 'inp', 'POSCAR')).structure
    parameters = IncarParser(file_path=data_path('test_converge_wc/pw/200', 'inp', 'INCAR')).incar
    parameters['system'] = 'test-case:test_converge_wc'
    parameters = {k: v for k, v in parameters.items() if k not in ['isif', 'ibrion', 'encut', 'nsw']}
    kpoints = KpointsParser(file_path=data_path('test_converge_wc/pw/200', 'inp', 'KPOINTS')).kpoints

    restart_clean_workdir = get_data_node('bool', False)
    restart_clean_workdir.store()

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.kpoints = kpoints
    inputs.parameters = get_data_node('dict', dict={'incar': parameters})
    inputs.potential_family = get_data_node('str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('dict', dict=POTCAR_MAP)
    inputs.options = get_data_node('dict',
                                   dict={
                                       'withmpi': False,
                                       'queue_name': 'None',
                                       'resources': {
                                           'num_machines': 1,
                                           'num_mpiprocs_per_machine': 1
                                       },
                                       'max_wallclock_seconds': 3600
                                   })
    inputs.max_iterations = get_data_node('int', 1)
    inputs.clean_workdir = get_data_node('bool', False)
    relax = AttributeDict()
    converge = AttributeDict()
    relax.perform = get_data_node('bool', False)
    converge.relax = get_data_node('bool', False)
    converge.testing = get_data_node('bool', True)
    converge.compress = get_data_node('bool', False)
    converge.displace = get_data_node('bool', False)
    converge.pwcutoff_start = get_data_node('float', 200)
    converge.pwcutoff_step = get_data_node('float', 50)
    converge.pwcutoff_samples = get_data_node('int', 3)
    inputs.relax = relax
    inputs.converge = converge
    inputs.verbose = get_data_node('bool', True)
    results, node = run.get_node(workchain, **inputs)
    assert node.exit_status == 0
    assert 'converge' in results
    converge = results['converge']
    assert 'data' in converge
    try:
        conv_data = converge['data']['pw_regular']
    except KeyError:
        pytest.fail('Did not find pw_regular in converge.data')
    conv_data_test = [
        [200.0, -10.77974998, 0.0, None, 0.5984],
        [250.0, -10.80762044, 0.0, None, 0.5912],
        [300.0, -10.82261992, 0.0, None, 0.5876],
    ]
    compare_converge_datum(conv_data, conv_data_test)

    assert 'pwcutoff_recommended' in converge
    try:
        _encut = converge['pwcutoff_recommended'].value
        np.testing.assert_equal(_encut, 300)
    except AttributeError:
        pytest.fail('pwcutoff_recommended does not have the expected format')


def test_converge_wc_kgrid(fresh_aiida_env, potentials, mock_vasp):
    """Test convergence workflow using mock code."""
    from aiida.orm import Code
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run

    workchain = WorkflowFactory('vasp.converge')

    mock_vasp.store()
    create_authinfo(computer=mock_vasp.computer).store()

    structure = PoscarParser(file_path=data_path('test_converge_wc/kgrid/8_8_8', 'inp', 'POSCAR')).structure
    parameters = IncarParser(file_path=data_path('test_converge_wc/kgrid/8_8_8', 'inp', 'INCAR')).incar
    parameters['system'] = 'test-case:test_converge_wc'
    pwcutoff = parameters['encut']
    parameters = {
        'incar': {k: v for k, v in parameters.items() if k not in ['isif', 'ibrion', 'encut', 'nsw']},
    }

    restart_clean_workdir = get_data_node('bool', False)
    restart_clean_workdir.store()

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.parameters = get_data_node('dict', dict=parameters)
    inputs.potential_family = get_data_node('str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('dict', dict=POTCAR_MAP)
    inputs.options = get_data_node('dict', dict={
        'queue_name': 'None',
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 1
        },
    })
    inputs.max_iterations = get_data_node('int', 1)
    inputs.clean_workdir = get_data_node('bool', False)
    relax = AttributeDict()
    relax.perform = get_data_node('bool', False)
    inputs.relax = relax
    converge = AttributeDict()
    converge.pwcutoff = get_data_node('float', pwcutoff)
    converge.k_dense = get_data_node('float', 0.21)  # 10x10x10 mesh
    converge.k_course = get_data_node('float', 0.27)  # 8x8x8 mesh
    converge.k_samples = get_data_node('int', 1)
    converge.testing = get_data_node('bool', True)
    inputs.converge = converge
    inputs.verbose = get_data_node('bool', True)
    results, node = run.get_node(workchain, **inputs)

    assert node.exit_status == 0
    assert 'converge' in results
    converge = results['converge']
    assert 'data' in converge
    conv_data = converge['data']
    try:
        conv_data = converge['data']['kpoints_regular']
    except KeyError:
        pytest.fail('Did not find kpoints_regular in converge.data')
    # this workflow should check 8x8x8 and 10x10x10 meshs
    np.testing.assert_allclose(conv_data[0][:4], [8.0, 8.0, 8.0, pwcutoff])
    np.testing.assert_allclose(conv_data[1][:4], [10.0, 10.0, 10.0, pwcutoff])

    assert 'kpoints_recommended' in converge
    try:
        _kpoints, _ = converge['kpoints_recommended'].get_kpoints_mesh()
        np.testing.assert_equal(_kpoints, [10.0, 10.0, 10.0])
    except AttributeError:
        pytest.fail('kpoints_recommended does not have the expected format')


def test_converge_wc_on_failed(fresh_aiida_env, vasp_params, potentials, mock_vasp):
    """Test convergence workflow using mock code."""
    from aiida.orm import Code
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run

    workchain = WorkflowFactory('vasp.converge')

    mock_vasp.store()
    create_authinfo(computer=mock_vasp.computer).store()

    structure = PoscarParser(file_path=data_path('test_converge_wc/pw/200', 'inp', 'POSCAR')).structure
    parameters = IncarParser(file_path=data_path('test_converge_wc/pw/200', 'inp', 'INCAR')).incar
    parameters['system'] = 'test-case:test_converge_wc'
    parameters = {k: v for k, v in parameters.items() if k not in ['isif', 'ibrion', 'encut', 'nsw']}
    kpoints = KpointsParser(file_path=data_path('test_converge_wc/pw/200', 'inp', 'KPOINTS')).kpoints

    restart_clean_workdir = get_data_node('bool', False)
    restart_clean_workdir.store()

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.kpoints = kpoints
    inputs.parameters = get_data_node('dict', dict={'incar': parameters})
    inputs.potential_family = get_data_node('str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('dict', dict=POTCAR_MAP)
    inputs.options = get_data_node('dict',
                                   dict={
                                       'withmpi': False,
                                       'queue_name': 'None',
                                       'resources': {
                                           'num_machines': 1,
                                           'num_mpiprocs_per_machine': 1
                                       },
                                       'max_wallclock_seconds': 3600
                                   })
    inputs.max_iterations = get_data_node('int', 1)
    inputs.clean_workdir = get_data_node('bool', False)
    relax = AttributeDict()
    converge = AttributeDict()
    relax.perform = get_data_node('bool', False)
    converge.relax = get_data_node('bool', False)
    converge.testing = get_data_node('bool', True)
    converge.compress = get_data_node('bool', False)
    converge.displace = get_data_node('bool', False)
    # intetionally make a called RelaxWorkChain fail by not proving output files at 'test_data/test_converge_wc/pw'
    converge.pwcutoff_start = get_data_node('float', 0.0)
    converge.pwcutoff_samples = get_data_node('int', 1)
    inputs.relax = relax
    inputs.converge = converge
    inputs.verbose = get_data_node('bool', True)
    results, node = run.get_node(workchain, **inputs)

    assert node.exit_status == 0
    assert 'converge' in results
    converge = results['converge']
    assert 'data' in converge
    assert 'pwcutoff_recommended' in converge
