"""
Test submitting a VaspWorkChain.

This does not seem to work, for `submit` the daemon will not pick up the workchain
and `run` just seems to get stuck after a while.
"""

# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name, import-outside-toplevel
from __future__ import print_function

import numpy as np
import pytest
from aiida import orm
from aiida.cmdline.utils.common import get_calcjob_report, get_workchain_report
from aiida.common.exceptions import InputValidationError
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import run
from aiida.plugins import WorkflowFactory
from aiida.plugins.factories import DataFactory


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_vasp_wc(run_vasp_process):
    """Test submitting only, not correctness, with mocked vasp code."""
    results, node = run_vasp_process(process_type='workchain')
    assert node.exit_status == 0
    assert 'retrieved' in results
    assert 'misc' in results
    assert 'remote_folder' in results
    misc = results['misc'].get_dict()
    assert np.amax(np.linalg.norm(misc['stress'], axis=1)) == pytest.approx(22.8499295)
    assert misc['total_energies']['energy_extrapolated'] == pytest.approx(-14.16209692)

    assert 'stdout' in node.called[0].inputs.monitors
    assert 'loop_time' in node.called[0].inputs.monitors


def si_structure():
    """
    Setup a silicon structure in a displaced FCC setting
    """
    structure_data = DataFactory('core.structure')
    alat = 3.9
    lattice = np.array([[0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]]) * alat
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


def setup_vasp_workchain(structure, incar, nkpts, potcar_family_name, potcar_mapping, code=None):
    """
    Setup the inputs for a VaspWorkChain.
    """

    inputs = AttributeDict()

    inputs.structure = structure
    inputs.parameters = orm.Dict(dict={'incar': incar})

    kpoints = orm.KpointsData()
    kpoints.set_kpoints_mesh((nkpts, nkpts, nkpts))
    inputs.kpoints = kpoints

    inputs.potential_family = orm.Str(potcar_family_name)
    inputs.potential_mapping = orm.Dict(dict=potcar_mapping)
    inputs.options = orm.Dict(
        dict={
            'withmpi': False,
            'queue_name': 'None',
            'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1},
            'max_wallclock_seconds': 3600,
        },
    )
    inputs.settings = orm.Dict(dict={'parser_settings': {'add_structure': True}})

    # If code is not passed, use the mock code
    if code is None:
        mock = orm.load_code('mock-vasp@localhost')
        inputs.code = mock
    else:
        inputs.code = code
    return inputs


def test_vasp_incar_case(aiida_profile):
    """Test catching inconsistent incar key cases"""
    workchain = WorkflowFactory('vasp.vasp')
    builder = workchain.get_builder()

    with pytest.warns(UserWarning, match="Key 'ALGO' converted to 'algo'"):
        builder.parameters = orm.Dict({'incar': {'ALGO': 21}})

    parameters = orm.Dict({'incar': {'ALGO': 21}})
    parameters.store()

    with pytest.raises(InputValidationError, match='Case inconsistency found in'):
        builder.parameters = parameters


def test_vasp_wc_nelm(aiida_profile, upload_potcar, potcar_family_name, potcar_mapping, mock_vasp_strict):
    """Test with mocked vasp code for handling electronic convergence issues"""

    workchain = WorkflowFactory('vasp.vasp')

    mock_vasp_strict.store()

    # create_authinfo(computer=mock_vasp_strict.computer, store=True)

    inputs = setup_vasp_workchain(si_structure(), INCAR_ELEC_CONV, 8, potcar_family_name, potcar_mapping)
    inputs.verbose = orm.Bool(True)
    results, node = run.get_node(workchain, **inputs)

    called_nodes = list(node.called)
    called_nodes.sort(key=lambda x: x.ctime)

    print(get_workchain_report(node, 'DEBUG'))
    for child in called_nodes:
        print(get_calcjob_report(child))

    child = called_nodes[0]
    print(child.base.repository.get_object_content('INCAR'))
    print(child.base.repository.get_object_content('POSCAR'))
    print(child.base.repository.get_object_content('KPOINTS'))
    print(child.outputs.retrieved.base.repository.get_object_content('vasp_output'))
    print(child.outputs.retrieved.base.repository.list_object_names())
    print(child.outputs.misc.get_dict())
    print(child.exit_status)

    child = called_nodes[1]
    print(child.base.repository.get_object_content('INCAR'))
    print(child.base.repository.get_object_content('POSCAR'))
    print(child.base.repository.get_object_content('KPOINTS'))
    print(child.outputs.retrieved.base.repository.get_object_content('vasp_output'))
    print(child.outputs.retrieved.base.repository.list_object_names())
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


@pytest.mark.parametrize(
    'incar,nkpts,exit_codes',
    [[INCAR_IONIC_CONV, 8, [702, 0]], [INCAR_IONIC_UNFINISHED, 16, [700, 0]]],
)
def test_vasp_wc_ionic_continue(
    aiida_profile, upload_potcar, potcar_family_name, potcar_mapping, mock_vasp_strict, incar, nkpts, exit_codes
):
    """Test with mocked vasp code for handling ionic convergence issues"""
    workchain = WorkflowFactory('vasp.vasp')

    mock_vasp_strict.store()
    # create_authinfo(computer=mock_vasp_strict.computer, store=True)

    inputs = setup_vasp_workchain(si_structure(), incar, nkpts, potcar_family_name, potcar_mapping)
    inputs.verbose = orm.Bool(True)
    # The test calculation contain NELM breaches during the relaxation - set to ignore it.
    inputs.handler_overrides = orm.Dict(dict={'ignore_nelm_breach_relax': True})
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


@pytest.mark.skip(reason='This test is not working yet')
def test_vasp_wc_ionic_magmom_carry(aiida_profile, upload_potcar, potcar_family_name, potcar_mapping, mock_vasp_strict):
    """Test with mocked vasp code for handling ionic convergence issues"""

    workchain = WorkflowFactory('vasp.vasp')

    mock_vasp_strict.store()
    # create_authinfo(computer=mock_vasp_strict.computer, store=True)

    incar = dict(INCAR_IONIC_CONV)
    incar['ispin'] = 2
    incar['lorbit'] = 10
    incar['nupdown'] = 2
    inputs = setup_vasp_workchain(si_structure(), incar, 8, potcar_family_name, potcar_mapping)
    inputs.verbose = orm.Bool(True)

    # The test calculation contain NELM breaches during the relaxation - set to ignore it.
    inputs.handler_overrides = orm.Dict(dict={'ignore_nelm_breach_relax': True})
    inputs.settings = orm.Dict(
        dict={
            'parser_settings': {
                'add_structure': True,
                'add_site_magnetization': True,
            }
        },
    )
    inputs.max_iterations = orm.Int(2)

    _, node = run.get_node(workchain, **inputs)
    assert node.exit_status == 0

    called_nodes = list(node.called)
    called_nodes.sort(key=lambda x: x.ctime)
    # Check that the second node takes the magnetization of the first node
    assert called_nodes[1].inputs.parameters['magmom'] == [0.646]
