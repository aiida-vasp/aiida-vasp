"""
Test submitting a RelaxWorkChain.

This does not seem to work, for `submit` the daemon will not pick up the workchain
and `run` just seems to get stuck after a while.
"""
# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name,no-member, import-outside-toplevel
from __future__ import print_function

import numpy as np
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida_vasp.parsers.content_parsers.incar import IncarParser
from aiida_vasp.parsers.content_parsers.kpoints import KpointsParser
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.parsers.node_composer import NodeComposer
from aiida_vasp.utils.aiida_utils import create_authinfo, get_data_node
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.mock_code import VaspMockRegistry


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
    global POTCAR_FAMILY_NAME  # noqa: PLW0603
    POTCAR_FAMILY_NAME = 'TEMP'
    potcar_data_cls = DataFactory('vasp.potcar')
    potcar_data_cls.upload_potcar_family(path, 'TEMP', 'TEMP-REALPOTCARS', stop_if_existing=False, dry_run=False)


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


def test_relax_wc(fresh_aiida_env, vasp_params, potentials, mock_vasp):
    # def test_relax_wc(fresh_aiida_env, vasp_params, potentials, mock_vasp, mock_relax_wc):
    """Test submitting only, not correctness, with mocked vasp code."""

    workchain = WorkflowFactory('vasp.relax')

    mock_vasp.store()
    create_authinfo(computer=mock_vasp.computer, store=True)

    structure = None
    with open(data_path('test_relax_wc', 'inp', 'POSCAR'), 'r', encoding='utf8') as handler:
        structure_parser = PoscarParser(handler=handler)
        structure = structure_parser.get_quantity('poscar-structure')
        structure = NodeComposer.compose_core_structure('core.structure', {'structure': structure})

    kpoints = None
    with open(data_path('test_relax_wc', 'inp', 'KPOINTS'), 'r', encoding='utf8') as handler:
        kpoints_parser = KpointsParser(handler=handler)
        kpoints = kpoints_parser.get_quantity('kpoints-kpoints')
        kpoints = NodeComposer.compose_core_array_kpoints('core.array.kpoints', {'kpoints': kpoints})
        kpoints.set_cell_from_structure(structure)

    parameters = None
    with open(data_path('test_relax_wc', 'inp', 'INCAR'), 'r', encoding='utf8') as handler:
        incar_parser = IncarParser(handler=handler)
        parameters = incar_parser.get_quantity('incar')

    parameters['system'] = 'test-case:test_relax_wc'
    parameters = {'incar': {k: v for k, v in parameters.items() if k not in ['isif', 'ibrion', 'nsw', 'ediffg']}}
    parameters['relax'] = {}
    parameters['relax']['perform'] = True
    parameters['relax']['algo'] = 'cg'
    parameters['relax']['force_cutoff'] = 0.01

    inputs = AttributeDict()
    inputs.code = orm.Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.kpoints = kpoints
    inputs.parameters = get_data_node('core.dict', dict=parameters)
    inputs.potential_family = get_data_node('core.str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('core.dict', dict=POTCAR_MAP)
    inputs.options = get_data_node(
        'core.dict',
        dict={
            'withmpi': False,
            'queue_name': 'None',
            'max_wallclock_seconds': 1,
            'import_sys_environment': True,
            'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1},
        },
    )
    inputs.max_iterations = get_data_node('core.int', 1)
    inputs.clean_workdir = get_data_node('core.bool', False)
    inputs.verbose = get_data_node('core.bool', True)
    results, node = run.get_node(workchain, **inputs)
    assert node.exit_status == 0
    assert 'relax' in results
    relax = results['relax']
    assert 'structure' in relax
    sites = relax['structure'].sites
    assert sites[0].kind_name == 'Si'
    assert sites[1].kind_name == 'Si'
    np.testing.assert_allclose(sites[0].position, [4.8125, 4.8125, 4.8125])
    np.testing.assert_allclose(sites[1].position, [0.6875, 0.6875, 0.715])


def test_mode_values():
    """Test that geometry relaxation modes either return a value or raise a ValueError"""
    from aiida_vasp.workchains.relax import RelaxWorkChain

    for pos in (True, False):
        for shape in (True, False):
            for volume in (True, False):
                try:
                    RelaxWorkChain.ModeEnum.get_from_dof(**{'positions': pos, 'shape': shape, 'volume': volume})
                except ValueError:
                    pass
                except Exception as exc:  # pylint: disable=broad-except
                    assert (
                        'Get from DOF function has to either return the correct value or raise a ValueError '
                        'for invalid combinations, instead got {} exception'.format(type(exc))
                    )


INCAR_SI_MAGMOM = {
    'encut': 240,
    'ismear': 0,
    'sigma': 0.1,
    'ediff': 1e-5,
    'nelm': 15,
    'potim': 0.1,
    'nupdown': 1,
}


def setup_vasp_relax_workchain(structure, incar, nkpts, code=None):
    """
    Setup the inputs for a VaspWorkChain.
    """
    # upload_real_pseudopotentials('/home/bonan/appdir/VASP/POTCARS/potpaw_PBE.54-2015_subset')
    inputs = AttributeDict()

    inputs.structure = structure
    inputs.parameters = get_data_node('core.dict', dict={'incar': incar})

    kpoints = get_data_node('core.array.kpoints')
    kpoints.set_kpoints_mesh((nkpts, nkpts, nkpts))
    inputs.kpoints = kpoints

    inputs.potential_family = get_data_node('core.str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('core.dict', dict=POTCAR_MAP)
    inputs.options = get_data_node(
        'core.dict',
        dict={
            'withmpi': False,
            'queue_name': 'None',
            'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1},
            'max_wallclock_seconds': 3600,
        },
    )
    inputs.settings = get_data_node('core.dict', dict={'parser_settings': {'add_structure': True}})
    inputs.relax = AttributeDict()
    inputs.relax.force_cutoff = orm.Float(5e-2)
    inputs.relax.volume = orm.Bool(True)
    inputs.relax.shape = orm.Bool(True)
    inputs.relax.steps = orm.Int(25)
    inputs.relax.perform = orm.Bool(True)
    inputs.relax.keep_magnetization = orm.Bool(True)

    # If code is not passed, use the mock code
    if code is None:
        mock = orm.Code.get_from_string('mock-vasp-strict@localhost')
        inputs.code = mock
    else:
        inputs.code = code
    return inputs


def test_vasp_wc_ionic_magmom_carry(fresh_aiida_env, potentials, mock_vasp_strict):
    """Test with mocked vasp code for handling ionic convergence issues"""

    workchain = WorkflowFactory('vasp.relax')

    mock_vasp_strict.store()
    create_authinfo(computer=mock_vasp_strict.computer, store=True)

    incar = dict(INCAR_SI_MAGMOM)
    incar['ispin'] = 2
    incar['lorbit'] = 10
    incar['nupdown'] = 2
    inputs = setup_vasp_relax_workchain(si_structure(), incar, 8)
    inputs.verbose = get_data_node('core.bool', True)

    # The test calculation contain NELM breaches during the relaxation - set to ignore it.
    inputs.settings = get_data_node(
        'core.dict',
        dict={
            'parser_settings': {
                'add_structure': True,
                'add_site_magnetization': True,
            }
        },
    )
    inputs.max_iterations = get_data_node('core.int', 3)

    _, node = run.get_node(workchain, **inputs)
    assert node.exit_status == 0

    called_nodes = list(node.called)
    called_nodes.sort(key=lambda x: x.ctime)
    # Check that the second node takes the magnetization of the first node
    assert called_nodes[1].inputs.site_magnetization['site_magnetization']['sphere']['x']['site_moment']['1']['tot'] == 0.643
    assert called_nodes[1].called[0].inputs.parameters['magmom'] == [0.643]
    assert 'site_magnetization' in node.outputs
    # upload_real_workchain(node, 'relax-wc-keep-magmom')
