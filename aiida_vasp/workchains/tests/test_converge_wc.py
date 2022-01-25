"""
Test submitting a ConvergenceWorkChain.

Only `run` currently works.

"""
# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name, too-many-statements, import-outside-toplevel, too-many-locals
from __future__ import print_function

import numpy as np
import pytest
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_node
from aiida_vasp.parsers.node_composer import NodeComposer
from aiida_vasp.parsers.content_parsers.kpoints import KpointsParser
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.parsers.content_parsers.incar import IncarParser
from aiida_vasp.utils.aiida_utils import create_authinfo


def test_converge_wc(fresh_aiida_env, potentials, mock_vasp):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import Code
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run

    workchain = WorkflowFactory('vasp.converge')

    mock_vasp.store()
    create_authinfo(computer=mock_vasp.computer, store=True)

    structure = None
    with open(data_path('test_converge_wc', 'inp', 'POSCAR'), 'r') as handler:
        structure_parser = PoscarParser(handler=handler)
        structure = structure_parser.get_quantity('poscar-structure')
        structure = NodeComposer.compose_structure('structure', {'structure': structure})
    parameters = None
    with open(data_path('test_converge_wc', 'inp', 'INCAR')) as handler:
        incar_parser = IncarParser(handler=handler)
        parameters = incar_parser.get_quantity('incar')

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
        conv_data.get_array('pw_regular')
    except KeyError:
        pytest.fail('Did not find pw_regular in converge.data')
    try:
        conv_data.get_array('kpoints_regular')
    except KeyError:
        pytest.fail('Did not find kpoints_regular in converge.data')

    assert 'pwcutoff_recommended' in converge
    try:
        _encut = converge['pwcutoff_recommended'].value
    except AttributeError:
        pytest.fail('pwcutoff_recommended does not have the expected format')
    assert 'kpoints_recommended' in converge
    try:
        _kpoints = converge['kpoints_recommended'].get_kpoints_mesh()
    except AttributeError:
        pytest.fail('kpoints_recommended does not have the expected format')


def test_converge_wc_pw(fresh_aiida_env, vasp_params, potentials, mock_vasp):
    """Test convergence workflow using mock code."""
    from aiida.orm import Code
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run

    workchain = WorkflowFactory('vasp.converge')

    mock_vasp.store()
    create_authinfo(computer=mock_vasp.computer).store()

    structure = None
    with open(data_path('test_converge_wc/pw/200', 'inp', 'POSCAR'), 'r') as handler:
        structure_parser = PoscarParser(handler=handler)
        structure = structure_parser.get_quantity('poscar-structure')
        structure = NodeComposer.compose_structure('structure', {'structure': structure})

    kpoints = None
    with open(data_path('test_converge_wc/pw/200', 'inp', 'KPOINTS'), 'r') as handler:
        kpoints_parser = KpointsParser(handler=handler)
        kpoints = kpoints_parser.get_quantity('kpoints-kpoints')
        kpoints = NodeComposer.compose_array_kpoints('array.kpoints', {'kpoints': kpoints})
        kpoints.set_cell_from_structure(structure)

    parameters = None
    with open(data_path('test_converge_wc/pw/200', 'inp', 'INCAR')) as handler:
        incar_parser = IncarParser(handler=handler)
        parameters = incar_parser.get_quantity('incar')

    parameters['system'] = 'test-case:test_converge_wc'
    parameters = {k: v for k, v in parameters.items() if k not in ['isif', 'ibrion', 'encut', 'nsw']}

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
    converge.pwcutoff_samples = get_data_node('int', 3)
    converge.k_samples = get_data_node('int', 3)
    inputs.relax = relax
    inputs.converge = converge
    inputs.verbose = get_data_node('bool', True)
    results, node = run.get_node(workchain, **inputs)
    assert node.exit_status == 0
    assert 'converge' in results
    converge = results['converge']
    assert 'data' in converge
    conv_data = converge['data']
    try:
        conv_data = conv_data.get_array('pw_regular')
    except KeyError:
        pytest.fail('Did not find pw_regular in converge.data')
    conv_data_test = np.array([[200.0, -10.77974998, 0.0, 0.0, 0.5984], [250.0, -10.80762044, 0.0, 0.0, 0.5912],
                               [300.0, -10.82261992, 0.0, 0.0, 0.5876]])
    np.testing.assert_allclose(conv_data, conv_data_test)

    assert 'pwcutoff_recommended' in converge
    try:
        _encut = converge['pwcutoff_recommended'].value
        np.testing.assert_equal(_encut, 300)
    except AttributeError:
        pytest.fail('pwcutoff_recommended does not have the expected format')
