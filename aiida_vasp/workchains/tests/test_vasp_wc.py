"""
Test submitting a VaspWorkChain.

This does not seem to work, for `submit` the daemon will not pick up the workchain
and `run` just seems to get stuck after a while.
"""
# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name, import-outside-toplevel
from __future__ import print_function

import pytest
import numpy as np
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.aiida_utils import get_data_node, aiida_version, cmp_version, create_authinfo


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_wc(fresh_aiida_env, vasp_params, potentials, vasp_kpoints, vasp_structure, mock_vasp):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import Code
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run

    workchain = WorkflowFactory('vasp.vasp')

    mock_vasp.store()
    create_authinfo(computer=mock_vasp.computer, store=True)

    kpoints, _ = vasp_kpoints
    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = vasp_structure
    inputs.parameters = get_data_node('dict', dict={'incar': vasp_params.get_dict()})
    inputs.kpoints = kpoints
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
    inputs.verbose = get_data_node('bool', True)
    results, node = run.get_node(workchain, **inputs)

    assert node.exit_status == 0
    assert 'retrieved' in results
    assert 'misc' in results
    assert 'remote_folder' in results
    misc = results['misc'].get_dict()
    assert misc['maximum_stress'] == 22.8499295
    assert misc['total_energies']['energy_extrapolated'] == -14.16209692


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_wc_chgcar(fresh_aiida_env, vasp_params, potentials, vasp_kpoints, vasp_structure, mock_vasp):
    """Test submitting only, not correctness, with mocked vasp code, test fetching of the CHGCAR."""
    from aiida.orm import Code
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run

    workchain = WorkflowFactory('vasp.vasp')

    mock_vasp.store()
    create_authinfo(computer=mock_vasp.computer, store=True)

    kpoints, _ = vasp_kpoints
    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = vasp_structure
    inputs.parameters = get_data_node('dict', dict={'incar': vasp_params.get_dict()})
    inputs.kpoints = kpoints
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
    inputs.settings = get_data_node('dict', dict={'ADDITIONAL_RETRIEVE_LIST': ['CHGCAR'], 'parser_settings': {'add_chgcar': True}})
    inputs.max_iterations = get_data_node('int', 1)
    inputs.clean_workdir = get_data_node('bool', False)
    inputs.verbose = get_data_node('bool', True)
    results, node = run.get_node(workchain, **inputs)
    assert node.exit_status == 0
    assert 'chgcar' in results
    assert results['chgcar'].get_content() == 'This is a test CHGCAR file.\n'


### COMPLEX WORHCAIN TEST ###


def si_structure():
    """
    Setup a silicon structure in a displaced FCC setting
    """
    from aiida.plugins import DataFactory

    structure_data = DataFactory('structure')
    alat = 3.9
    lattice = np.array([[.5, .5, 0], [0, .5, .5], [.5, 0, .5]]) * alat
    structure = structure_data(cell=lattice)
    positions = [[0.1, 0.0, 0.0]]
    for pos_direct in positions:
        pos_cartesian = np.dot(pos_direct, lattice)
        structure.append_atom(position=pos_cartesian, symbols='Si')
    return structure


# TEST INPUT FOR AUTOMATIC correction of NELM
# calculation should finish in the second run where the calculation
INCAR_ELEC_CONV = {
    'encut': 240,
    'ismear': 0,
    'sigma': 0.1,
    'ediff': 1e-9,
    'nelm': 7,
    'ibrion': -1,
    'potim': 0.01,
    'nsw': -1,
    'isif': 3,
    # 'ediffg': -0.01
}

INCAR_IONIC_CONV = {
    'encut': 240,
    'ismear': 0,
    'sigma': 0.1,
    'ediff': 1e-9,
    'nelm': 15,
    'ibrion': 1,
    'potim': 0.1,
    'nsw': 5,
    'isif': 3,
}

# Parameters for test handling unfinished VASP. The first iteration was killed manually.
INCAR_IONIC_UNFINISHED = {
    'encut': 500,
    'ismear': 0,
    'isym': 0,
    'sigma': 0.1,
    'ediff': 1e-9,
    'nelm': 15,
    'ibrion': 1,
    'potim': 0.1,
    'nsw': 20,
    'isif': 3,
}


def setup_vasp_workchain(structure, incar, nkpts):
    """
    Setup the inputs for a VaspWorkChain.
    """
    from aiida.orm import Code

    inputs = AttributeDict()

    inputs.structure = structure
    inputs.parameters = get_data_node('dict', dict={'incar': incar})

    kpoints = get_data_node('array.kpoints')
    kpoints.set_kpoints_mesh((nkpts, nkpts, nkpts))
    inputs.kpoints = kpoints

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
    inputs.settings = get_data_node('dict', dict={'parser_settings': {'add_structure': True}})

    mock = Code.get_from_string('mock-vasp-strict@localhost')
    inputs.code = mock
    return inputs


def test_vasp_wc_nelm(fresh_aiida_env, potentials, mock_vasp_strict):
    """Test with mocked vasp code for handling electronic convergence issues"""
    from aiida.orm import Code
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run
    from aiida.cmdline.utils.common import get_calcjob_report, get_workchain_report

    workchain = WorkflowFactory('vasp.vasp')

    mock_vasp_strict.store()
    create_authinfo(computer=mock_vasp_strict.computer, store=True)

    inputs = setup_vasp_workchain(si_structure(), INCAR_ELEC_CONV, 8)
    inputs.verbose = get_data_node('bool', True)
    results, node = run.get_node(workchain, **inputs)

    called_nodes = list(node.called)
    called_nodes.sort(key=lambda x: x.ctime)

    print(get_workchain_report(node, 'DEBUG'))
    for child in called_nodes:
        print(get_calcjob_report(child))

    child = called_nodes[0]
    print(child.get_object_content('INCAR'))
    print(child.get_object_content('POSCAR'))
    print(child.get_object_content('KPOINTS'))
    print(child.outputs.retrieved.get_object_content('vasp_output'))
    print(child.outputs.retrieved.list_object_names())
    print(child.outputs.misc.get_dict())
    print(child.exit_status)

    child = called_nodes[1]
    print(child.get_object_content('INCAR'))
    print(child.get_object_content('POSCAR'))
    print(child.get_object_content('KPOINTS'))
    print(child.outputs.retrieved.get_object_content('vasp_output'))
    print(child.outputs.retrieved.list_object_names())
    print(child.outputs.misc.get_dict())
    print(child.exit_status)

    assert node.exit_status == 0
    assert 'retrieved' in results
    assert 'misc' in results
    assert 'remote_folder' in results

    assert results['misc']['total_energies']['energy_extrapolated'] == -4.82467802

    # Sort the called nodes by creation time
    called_nodes = list(node.called)
    called_nodes.sort(key=lambda x: x.ctime)

    assert called_nodes[0].exit_status == 701
    assert called_nodes[1].exit_status == 0


@pytest.mark.parametrize('incar,nkpts,exit_codes', [[INCAR_IONIC_CONV, 8, [702, 0]], [INCAR_IONIC_UNFINISHED, 16, [700, 0]]])
def test_vasp_wc_ionic_continue(fresh_aiida_env, potentials, mock_vasp_strict, incar, nkpts, exit_codes):
    """Test with mocked vasp code for handling ionic convergence issues"""
    from aiida.orm import Code
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run

    workchain = WorkflowFactory('vasp.vasp')

    mock_vasp_strict.store()
    create_authinfo(computer=mock_vasp_strict.computer, store=True)

    inputs = setup_vasp_workchain(si_structure(), incar, nkpts)
    inputs.verbose = get_data_node('bool', True)
    results, node = run.get_node(workchain, **inputs)

    assert node.exit_status == 0
    assert 'retrieved' in results
    assert 'misc' in results
    assert 'remote_folder' in results

    assert results['misc']['run_status']['ionic_converged']

    # Sort the called nodes by creation time
    called_nodes = list(node.called)
    called_nodes.sort(key=lambda x: x.ctime)

    # Check the child status - here the first calculation is not finished but the second one is
    for idx, code in enumerate(exit_codes):
        assert called_nodes[idx].exit_status == code
