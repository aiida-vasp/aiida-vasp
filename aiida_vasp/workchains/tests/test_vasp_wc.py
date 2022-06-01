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
from aiida.manage.tests.pytest_fixtures import aiida_caplog
from aiida.plugins.factories import DataFactory

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.aiida_utils import get_data_node, aiida_version, cmp_version, create_authinfo
from aiida_vasp.utils.mock_code import VaspMockRegistry


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_wc(fresh_aiida_env, run_vasp_process):
    """Test submitting only, not correctness, with mocked vasp code."""
    results, node = run_vasp_process(process_type='workchain')
    assert node.exit_status == 0
    assert 'retrieved' in results
    assert 'misc' in results
    assert 'remote_folder' in results
    misc = results['misc'].get_dict()
    assert misc['maximum_stress'] == pytest.approx(22.8499295)
    assert misc['total_energies']['energy_extrapolated'] == pytest.approx(-14.16209692)


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_wc_chgcar(fresh_aiida_env, run_vasp_process, aiida_caplog):
    """Test submitting only, not correctness, with mocked vasp code, test fetching and parsing of the CHGCAR content."""
    settings = {'ADDITIONAL_RETRIEVE_LIST': ['CHGCAR'], 'parser_settings': {'add_charge_density': True}}
    results, node = run_vasp_process(settings=settings, process_type='workchain')
    assert node.exit_status == 0
    assert 'charge_density' in results
    assert 'misc' in results
    test_array = np.array([[[0.09329446, 0.18658892, 0.27988338], [0.37317784, 0.4664723, 0.55976676], [0.65306122, 0.74635569, 0.83965015],
                            [0.93294461, 1.02623907, 1.11953353]],
                           [[1.21282799, 1.30612245, 1.39941691], [1.49271137, 1.58600583, 1.67930029],
                            [1.77259475, 1.86588921, 1.95918367], [2.05247813, 2.14577259, 2.23906706]],
                           [[2.33236152, 2.42565598, 2.51895044], [2.6122449, 2.70553936, 2.79883382], [2.89212828, 2.98542274, 3.0787172],
                            [3.17201166, 3.26530612, 3.35860058]],
                           [[3.45189504, 3.5451895, 3.63848397], [3.73177843, 3.82507289, 3.91836735], [4.01166181, 4.10495627, 4.19825073],
                            [4.29154519, 4.38483965, 4.47813411]],
                           [[4.57142857, 4.66472303, 4.75801749], [4.85131195, 4.94460641, 5.03790087], [5.13119534, 5.2244898, 5.31778426],
                            [5.41107872, 5.50437318, 5.59766764]]])
    charge_density = results['charge_density'].get_array('charge_density')
    assert np.allclose(charge_density, test_array)


def upload_real_workchain(node, name):
    """
    Upload the workchain to the repository to make it work with mocking

    This function should be called once after the REAL vasp calculation is run during the test
    """
    reg = VaspMockRegistry()
    print(reg.base_path)
    reg.upload_aiida_work(node, name)


def upload_real_pseudopotentials(path):
    """
    Upload real pseudopotentials for workchain test mock deposition


    This function should be called once before the REAL vasp calculation is launch to setup the
    correct POTCARs
    """
    global POTCAR_FAMILY_NAME  # pylint: disable=global-statement
    POTCAR_FAMILY_NAME = 'TEMP'
    potcar_data_cls = DataFactory('vasp.potcar')
    potcar_data_cls.upload_potcar_family(path, 'TEMP', 'TEMP-REALPOTCARS', stop_if_existing=False, dry_run=False)


### COMPLEX WORKCHAIN TEST ###


def si_structure():
    """
    Setup a silicon structure in a displaced FCC setting
    """
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


def setup_vasp_workchain(structure, incar, nkpts, code=None):
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

    # If code is not passed, use the mock code
    if code is None:
        mock = Code.get_from_string('mock-vasp-strict@localhost')
        inputs.code = mock
    else:
        inputs.code = code
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

    assert results['misc']['total_energies']['energy_extrapolated'] == pytest.approx(-4.82467802)

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
    # The test calculation contain NELM breaches during the relaxation - set to ignore it.
    inputs.handler_overrides = get_data_node('dict', dict={'ignore_nelm_breach_relax': True})
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


def test_vasp_wc_ionic_magmom_carry(fresh_aiida_env, potentials, mock_vasp_strict):
    """Test with mocked vasp code for handling ionic convergence issues"""
    from aiida.orm import Code
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run

    workchain = WorkflowFactory('vasp.vasp')

    mock_vasp_strict.store()
    create_authinfo(computer=mock_vasp_strict.computer, store=True)

    incar = dict(INCAR_IONIC_CONV)
    incar['ispin'] = 2
    incar['lorbit'] = 10
    incar['nupdown'] = 2
    inputs = setup_vasp_workchain(si_structure(), incar, 8)
    inputs.verbose = get_data_node('bool', True)

    # The test calculation contain NELM breaches during the relaxation - set to ignore it.
    inputs.handler_overrides = get_data_node('dict', dict={'ignore_nelm_breach_relax': True})
    inputs.settings = get_data_node('dict', dict={'parser_settings': {
        'add_structure': True,
        'add_site_magnetization': True,
    }})
    inputs.max_iterations = get_data_node('int', 2)

    _, node = run.get_node(workchain, **inputs)
    assert node.exit_status == 0

    called_nodes = list(node.called)
    called_nodes.sort(key=lambda x: x.ctime)
    # Check that the second node takes the magnetization of the first node
    assert called_nodes[1].inputs.parameters['magmom'] == [0.646]
