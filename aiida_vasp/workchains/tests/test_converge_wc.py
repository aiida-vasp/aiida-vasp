"""
Test submitting a ConvergenceWorkChain.

Only `run` currently works.

"""
# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name, too-many-statements
from __future__ import print_function

import numpy as np
import pytest
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_node, get_current_user, aiida_version, cmp_version
from aiida_vasp.parsers.file_parsers.kpoints import KpointsParser
from aiida_vasp.parsers.file_parsers.poscar import PoscarParser
from aiida_vasp.parsers.file_parsers.incar import IncarParser
from aiida_vasp.utils.aiida_utils import create_authinfo


@pytest.mark.wc
@pytest.mark.skip(reason='We can not run consecutive tests of this type due to a bug in AiiDA '
                  'which will be fixed when the django_jsonb branch is merged. Disabling the '
                  'simplest test.')
# @pytest.mark.skip(reason='Travis has problems with runs that continues after 10 min. In addition, '
#                   'we cannnot run two executive tests, e.g. test_converge_wc and test_converge_pw '
#                   'as we get IntegrityError from Django. The source of this is currently unknown, '
#                   'so we disable the least detailed test.')
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
    parameters = {k: v for k, v in parameters.items() if k not in ['isif', 'ibrion', 'encut', 'nsw']}
    parameters['system'] = 'test-case:test_converge_wc'

    restart_clean_workdir = get_data_node('bool', False)
    restart_clean_workdir.store()

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.parameters = get_data_node('dict', dict=parameters)
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
    converge.converge_relax = get_data_node('bool', True)
    converge.compress = get_data_node('bool', False)
    converge.displace = get_data_node('bool', False)
    converge.encut_samples = get_data_node('int', 3)
    converge.k_samples = get_data_node('int', 3)
    relax.perform = get_data_node('bool', True)
    inputs.relax = relax
    inputs.converge = converge
    inputs.verbose = get_data_node('bool', True)
    results, node = run.get_node(workchain, **inputs)
    assert node.exit_status == 0
    assert 'convergence.data' in results
    assert 'structure_relaxed' in results

    conv_data = results['convergence.data']
    try:
        conv_data.get_array('pw_regular')
    except KeyError:
        pytest.fail('Did not find pw_regular in convergence.data')
    try:
        conv_data.get_array('kpoints_regular')
    except KeyError:
        pytest.fail('Did not find kpoints_regular in convergence.data')


@pytest.mark.skip(reason='Too slow')
@pytest.mark.wc
def test_converge_wc_pw(fresh_aiida_env, vasp_params, potentials, mock_vasp):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import Code
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run

    workchain = WorkflowFactory('vasp.converge')

    mock_vasp.store()
    create_authinfo(computer=mock_vasp.computer).store()

    structure = PoscarParser(file_path=data_path('test_converge_wc/pw/200', 'inp', 'POSCAR')).structure
    parameters = IncarParser(file_path=data_path('test_converge_wc/pw/200', 'inp', 'INCAR')).incar
    parameters = {k: v for k, v in parameters.items() if k not in ['isif', 'ibrion', 'encut', 'nsw']}
    kpoints = KpointsParser(file_path=data_path('test_converge_wc/pw/200', 'inp', 'KPOINTS')).kpoints
    parameters['system'] = 'test-case:test_converge_wc'

    restart_clean_workdir = get_data_node('bool', False)
    restart_clean_workdir.store()

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.kpoints = kpoints
    inputs.parameters = get_data_node('dict', dict=parameters)
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
    converge.converge_relax = get_data_node('bool', False)
    converge.testing = get_data_node('bool', True)
    converge.compress = get_data_node('bool', False)
    converge.displace = get_data_node('bool', False)
    converge.encut_samples = get_data_node('int', 3)
    converge.k_samples = get_data_node('int', 3)
    inputs.relax = relax
    inputs.converge = converge
    inputs.verbose = get_data_node('bool', True)
    results, node = run.get_node(workchain, **inputs)
    assert node.exit_status == 0
    assert 'convergence.data' in results
    conv_data = results['convergence.data']
    try:
        conv_data = conv_data.get_array('pw_regular')
    except KeyError:
        pytest.fail('Did not find pw_regular in convergence.data')
    conv_data_test = np.array([[200.0, -10.77974998, 0.0, 0.0, 0.5984], [250.0, -10.80762044, 0.0, 0.0, 0.5912],
                               [300.0, -10.82261992, 0.0, 0.0, 0.5876]])
    np.testing.assert_allclose(conv_data, conv_data_test)
