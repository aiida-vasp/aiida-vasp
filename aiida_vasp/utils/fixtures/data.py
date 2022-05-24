"""
Fixtures related to data.

-------------------------
Here different pytest fixtures are set up. They typically contain computers,
VASP inputs etc. which you would need to mock a VASP job in this plugin. It
also contains the set up of the mock VASP executable, which is used to test the
workchains.
"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name,too-many-function-args,
# pylint: disable=protected-access,abstract-class-instantiated,no-value-for-parameter,unexpected-keyword-arg, import-outside-toplevel
import os
from collections import OrderedDict
import subprocess as sp
from pathlib import Path

import numpy
import pytest
from pymatgen.io.vasp import Poscar

from aiida.orm import Computer, FolderData
from aiida.common.exceptions import NotExistent
from aiida.common.extendeddicts import AttributeDict
from aiida.manage.tests import TemporaryProfileManager
from aiida_vasp.utils.aiida_utils import get_data_node, get_data_class
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.parsers.content_parsers.incar import IncarParser
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.parsers.content_parsers.vasprun import VasprunParser
from aiida_vasp.parsers.content_parsers.outcar import OutcarParser, VtstNebOutcarParser
from aiida_vasp.parsers.content_parsers.doscar import DoscarParser
from aiida_vasp.parsers.content_parsers.eigenval import EigenvalParser
from aiida_vasp.parsers.content_parsers.chgcar import ChgcarParser
from aiida_vasp.parsers.content_parsers.kpoints import KpointsParser
from aiida_vasp.utils.general import copytree
from aiida_vasp.parsers.content_parsers.stream import StreamParser
from aiida_vasp.data.potcar import OLD_POTCAR_FAMILY_TYPE, PotcarGroup, Group

POTCAR_FAMILY_NAME = 'test_family'
POTCAR_MAP = {'In': 'In_sv', 'In_d': 'In_d', 'As': 'As', 'Ga': 'Ga', 'Si': 'Si', 'P': 'P', 'S': 'S', 'Zn': 'Zn', 'N': 'N', 'H': 'H'}


def path_file_and_settings(name, param):
    """Locate folder, filename and settings from param. Return the path and settings."""
    settings = {}
    if isinstance(param, list):
        if len(param) == 3:
            folder, name, settings = param
        elif len(param) == 2:
            folder, name = param
        else:
            raise IndexError('Please supply either folder and name, or folder, name and settings to the parser fixtures')
    else:
        folder = param
    path = data_path(folder, name)

    return path, settings


@pytest.fixture(scope='session')
def localhost_dir(tmp_path_factory):
    return tmp_path_factory.mktemp('localhost_work')


@pytest.fixture
def localhost(fresh_aiida_env, localhost_dir):
    """Fixture for a local computer called localhost. This is currently not in the AiiDA fixtures."""
    try:
        computer = Computer.objects.get(label='localhost')
    except NotExistent:
        computer = Computer(label='localhost',
                            hostname='localhost',
                            transport_type='local',
                            scheduler_type='direct',
                            workdir=str(localhost_dir)).store()
        computer.set_minimum_job_poll_interval(0.)
    return computer


@pytest.fixture
def vasp_params(fresh_aiida_env):
    incar_data = get_data_class('dict')(dict={'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'sigma': 0.5, 'magmom': '30 * 2*0.'})
    return incar_data


@pytest.fixture
def vasp2w90_params(fresh_aiida_env, vasp_params):
    vasp_params_data = vasp_params()
    incar_data = get_data_class('dict')(dict=vasp_params_data.code.get_dict().update({'lwannier90': True}))
    return incar_data


@pytest.fixture
def potcar_node_pair(fresh_aiida_env):
    """Create a POTCAR node pair."""
    potcar_path = data_path('potcar', 'As', 'POTCAR')
    potcar_file_node = get_data_node('vasp.potcar_file', file=potcar_path)
    potcar_file_node.store()
    return {'file': potcar_file_node, 'potcar': get_data_class('vasp.potcar').find_one(symbol='As')}


@pytest.fixture
def temp_pot_folder(tmp_path):
    """A temporary copy of the potcar test data folder, to avoid extracting tar objects inside the repo."""
    potcar_ga = Path(data_path('potcar')) / 'Ga'
    assert not potcar_ga.exists()
    pot_archive = Path(data_path('potcar'))
    target = tmp_path / 'potentials'
    # Ensure that the target path exists
    Path(target).mkdir(exist_ok=True)
    copytree(pot_archive, target)
    return target


# pylint: disable=protected-access
def duplicate_potcar_data(potcar_node):
    """Create and store (and return) a duplicate of a given PotcarData node."""
    from aiida_vasp.data.potcar import temp_potcar
    file_node = get_data_node('vasp.potcar_file')
    with temp_potcar(potcar_node.get_content()) as potcar_file:
        file_node.add_file(potcar_file)
        file_node.set_attribute('sha512', 'abcd')
        file_node.set_attribute('full_name', potcar_node.full_name)
        file_node.store()
    data_node, _ = get_data_class('vasp.potcar').get_or_create(file_node)
    return data_node


@pytest.fixture
def potcar_family(fresh_aiida_env, temp_pot_folder):
    """Create a POTCAR family."""
    potcar_ga = Path(data_path('potcar')) / 'Ga'
    family_name = POTCAR_FAMILY_NAME
    family_desc = 'A POTCAR family used as a test fixture. Contains only unusable POTCAR files.'
    potcar_cls = get_data_class('vasp.potcar')
    potcar_cls.upload_potcar_family(str(temp_pot_folder), family_name, family_desc, stop_if_existing=False)
    if len(potcar_cls.find(full_name='In_d')) == 1:
        family_group = potcar_cls.get_potcar_group(POTCAR_FAMILY_NAME)
        in_d = potcar_cls.find(full_name='In_d')[0]
        in_d_double = duplicate_potcar_data(in_d)
        family_group.add_nodes(in_d_double)
        assert in_d.uuid == potcar_cls.find(full_name='In_d')[0].uuid
    assert 'As' in potcar_cls.get_full_names(POTCAR_FAMILY_NAME, 'As')
    assert 'Ga' in potcar_cls.get_full_names(POTCAR_FAMILY_NAME, 'Ga')
    assert 'In_d' in potcar_cls.get_full_names(POTCAR_FAMILY_NAME, 'In')
    assert not potcar_ga.exists()
    return family_name


@pytest.fixture
def legacy_potcar_family(potcar_family):
    """
    Fixture from creating an legacy potcar group

    Returns a tuple of group label and the LegacyPotcarGroup with the old type_string
    """

    class LegacyPotcarGroup(Group):
        """Old style group with the old type string"""

    # Override the _type_string class property which is supposed to be loaded from the entrypoint.
    LegacyPotcarGroup._type_string = OLD_POTCAR_FAMILY_TYPE
    new_group = PotcarGroup.objects.get(label=potcar_family)
    old_group = LegacyPotcarGroup(label=potcar_family + '_migrate_test')
    old_group.store()

    # Add the nodes from the new group to the old group
    old_group.add_nodes(list(new_group.nodes))
    return old_group.label, LegacyPotcarGroup


@pytest.fixture
def potentials(potcar_family):
    """Fixture for two incomplete POTPAW potentials."""
    potcar_cls = get_data_class('vasp.potcar')
    potentials = potcar_cls.get_potcars_dict(['In', 'In_d', 'As'], family_name=potcar_family, mapping=POTCAR_MAP)

    return potentials


@pytest.fixture(params=['cif', 'str'])
def vasp_structure(request, fresh_aiida_env):
    """Fixture: StructureData or CifData."""
    from aiida.plugins import DataFactory
    if request.param == 'cif':
        cif_path = data_path('cif', 'EntryWithCollCode43360.cif')
        structure = DataFactory('cif').get_or_create(cif_path)[0]
    elif request.param == 'str':
        larray = numpy.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
        alat = 6.058
        structure = DataFactory('structure')(cell=larray * alat)
        structure.append_atom(position=[0, 0, 0], symbols='In')
        structure.append_atom(position=[.25, .25, .25], symbols='As')
        structure.append_atom(position=[.25, .33, .34], symbols='As')
        structure.append_atom(position=[.5, .5, .5], symbols='In', name='In_d')
        structure.append_atom(position=[.7896, .6234, .5], symbols='In', name='In_d')
        structure.append_atom(position=[.75, .75, .75], symbols='As')
    elif request.param == 'str-Al':
        larray = numpy.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        alat = 4.04
        structure = DataFactory('structure')(cell=larray * alat)
        structure.append_atom(position=numpy.array([0, 0, 0]) * alat, symbols='Al')
        structure.append_atom(position=numpy.array([0, .5, .5]) * alat, symbols='Al')
        structure.append_atom(position=numpy.array([.5, 0, .5]) * alat, symbols='Al')
        structure.append_atom(position=numpy.array([.5, .5, 0]) * alat, symbols='Al')
    elif request.param == 'str-InAs':
        structure_cls = DataFactory('structure')
        structure = structure_cls(cell=numpy.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * 6.058)
        structure.append_atom(position=(0, 0, 0), symbols='In', name='Hamburger')
        structure.append_atom(position=(0.25, 0.25, 0.25), symbols='As', name='Pizza')
    return structure


@pytest.fixture()
def vasp_structure_poscar(vasp_structure):
    """Fixture: Well formed POSCAR contents."""
    aiida_structure = vasp_structure
    if isinstance(vasp_structure, get_data_class('cif')):
        ase_structure = vasp_structure.get_ase()
        aiida_structure = get_data_node('structure', ase=ase_structure)
    writer = PoscarParser(data=aiida_structure)
    return writer


@pytest.fixture(params=['mesh', 'list'])
def vasp_kpoints(request, fresh_aiida_env):
    """Fixture: (kpoints object, resulting KPOINTS)."""
    from aiida.plugins import DataFactory
    if request.param == 'mesh':
        kpoints = DataFactory('array.kpoints')()
        kpoints.set_kpoints_mesh([2, 2, 2])
        ref_kpoints = _ref_kp_mesh()
    elif request.param == 'list':
        kpoints = DataFactory('array.kpoints')()
        kpoints.set_kpoints([[0., 0., 0.], [0., 0., .5]], weights=[1., 1.])
        ref_kpoints = _ref_kp_list()
    return kpoints, ref_kpoints


@pytest.fixture()
def vasp_inputs(fresh_aiida_env, vasp_params, vasp_kpoints, vasp_structure, potentials, vasp_code):
    """Inputs dictionary for CalcJob Processes."""
    from aiida.orm import Dict

    def inner(settings=None, parameters=None):

        inputs = AttributeDict()

        metadata = AttributeDict({'options': {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}})

        if settings is not None:
            inputs.settings = Dict(dict=settings)

        if isinstance(parameters, dict):
            parameters = get_data_class('dict')(dict=parameters)

        if parameters is None:
            parameters = AttributeDict(vasp_params.get_dict())
            parameters = get_data_class('dict')(dict=parameters)
        inputs.code = vasp_code
        inputs.metadata = metadata
        inputs.parameters = parameters
        inputs.kpoints, _ = vasp_kpoints
        inputs.structure = vasp_structure
        inputs.potential = potentials

        return inputs

    return inner


@pytest.fixture()
def vasp2w90_inputs(
        fresh_aiida_env,  # yapf: disable
        vasp_params,  # yapf: disable
        vasp_kpoints,  # yapf: disable
        vasp_structure,  # yapf: disable
        potentials,  # yapf: disable
        vasp_code,  # yapf: disable
        wannier_projections,  # yapf: disable
        wannier_params  # yapf: disable
):  # pylint: disable=too-many-arguments
    """Inputs dictionary for CalcJob Processes."""
    from aiida.orm import Dict

    def inner(settings=None, parameters=None):

        inputs = AttributeDict()

        metadata = AttributeDict({'options': {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}})

        if settings is not None:
            inputs.settings = Dict(dict=settings)

        if isinstance(parameters, dict):
            parameters = get_data_class('dict')(dict=parameters)

        if parameters is None:
            parameters = AttributeDict(vasp_params.get_dict())
            parameters = get_data_class('dict')(dict=parameters)

        inputs.code = vasp_code
        inputs.metadata = metadata
        inputs.parameters = parameters
        inputs.kpoints, _ = vasp_kpoints
        inputs.structure = vasp_structure
        inputs.potential = potentials

        inputs.wannier_parameters = wannier_params
        inputs.wannier_projections = wannier_projections

        return inputs

    return inner


@pytest.fixture()
def vasp_neb_inputs(fresh_aiida_env, vasp_params, vasp_kpoints, vasp_structure, potentials, vasp_code):
    """Inputs dictionary for CalcJob Processes."""
    from aiida.orm import Dict

    def inner(settings=None, parameters=None):

        inputs = AttributeDict()

        metadata = AttributeDict({'options': {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}})

        if settings is not None:
            inputs.settings = Dict(dict=settings)

        if isinstance(parameters, dict):
            parameters = get_data_class('dict')(dict=parameters)

        if parameters is None:
            parameters = AttributeDict(vasp_params.get_dict())
            parameters['images'] = 3
            parameters = get_data_class('dict')(dict=parameters)

        inputs.code = vasp_code
        inputs.metadata = metadata
        inputs.parameters = parameters
        inputs.kpoints, _ = vasp_kpoints

        inputs.initial_structure = vasp_structure
        inputs.final_structure = vasp_structure

        inputs.potential = potentials

        neb_images = {f'images_{idx:02d}': vasp_structure for idx in range(1, 4)}
        inputs.neb_images = neb_images

        return inputs

    return inner


@pytest.fixture()
def vasp_code(localhost):
    """Fixture for a vasp code, the executable it points to does not exist."""
    from aiida.orm import Code
    if not localhost.pk:
        localhost.store()
    code = Code()
    code.label = 'vasp'
    code.description = 'VASP code'
    code.set_remote_computer_exec((localhost, '/usr/local/bin/vasp'))
    code.set_input_plugin_name('vasp.vasp')
    return code


@pytest.fixture()
def mock_vasp(fresh_aiida_env, localhost):
    """
    Give an mock-up of the VASP executable

    This code will always create the output object even if no matching
    calculations from the registry is found. This makes it suitable for simple
    tests.
    """
    return _mock_vasp(fresh_aiida_env, localhost, 'mock-vasp')


@pytest.fixture()
def mock_vasp_strict(fresh_aiida_env, localhost):
    """
    Give an mock-up of the VASP executable with strict input matching.

    This code will not create the output object unless matching calculations from the
    registry is found. It is suitable for testsing complex multi-step workchains.
    tests.
    """
    return _mock_vasp(fresh_aiida_env, localhost, 'mock-vasp-strict')


def _mock_vasp(fresh_aiida_env, localhost, exec_name):
    """
    Points to a mock-up of a VASP executable.

    If environmental variable REAL_VASP_PATH is set, switch the code
    to point to the REAL VASP executable. This is used to generate the
    actual outputs for mock tests later
    """
    from aiida.orm import Code
    from aiida.orm.querybuilder import QueryBuilder
    query_builder = QueryBuilder()
    query_builder.append(Code, tag='code')
    query_builder.add_filter('code', {'label': {'==': exec_name}})
    query_results = query_builder.all()
    if query_results:
        code = query_results[0][0]
    else:
        os_env = os.environ.copy()
        if not localhost.pk:
            localhost.store()
        # returns unicode
        mock_vasp_path = sp.check_output(['which', exec_name], env=os_env, universal_newlines=True).strip()

        # Allow overriding mock using REAL code, this is used for running the actual
        # calculation once and deposit the results in the registry
        if os.environ.get('REAL_VASP_PATH'):
            mock_vasp_path = os.environ['REAL_VASP_PATH']

        code = Code()
        code.label = exec_name
        code.description = 'Mock VASP for tests'
        code.set_remote_computer_exec((localhost, mock_vasp_path))
        code.set_input_plugin_name('vasp.vasp')
        if isinstance(fresh_aiida_env._manager, TemporaryProfileManager):
            aiidapath = Path(fresh_aiida_env._manager.root_dir) / '.aiida'
        else:
            aiidapath = Path(os.environ.get('AIIDA_PATH', os.environ.get('HOME'))) / '/.aiida'
        code.set_prepend_text('export AIIDA_PATH={}'.format(aiidapath))
        code.store()
        code.set_extra('is_mock_code', True)

    return code


@pytest.fixture()
def vasp_chgcar(fresh_aiida_env):
    """CHGCAR node and reference fixture."""
    from aiida.plugins import DataFactory
    chgcar_path = data_path('chgcar', 'CHGCAR')
    chgcar = DataFactory('vasp.chargedensity')(file=chgcar_path)
    with open(chgcar_path, 'r') as ref_chgcar_fo:
        ref_chgcar = ref_chgcar_fo.read()
    return chgcar, ref_chgcar


@pytest.fixture()
def vasp_wavecar(fresh_aiida_env):
    """WAVECAR node and reference fixture."""
    from aiida.plugins import DataFactory
    wavecar_path = data_path('wavecar', 'WAVECAR')
    wavecar = DataFactory('vasp.wavefun')(file=wavecar_path)
    with open(wavecar_path, 'r') as ref_wavecar_fo:
        ref_wavecar = ref_wavecar_fo.read()
    return wavecar, ref_wavecar


@pytest.fixture
def ref_incar():
    with open(data_path('incar', 'INCAR'), 'r') as reference_incar_fo:
        # yield reference_incar_fo.read().strip()
        yield reference_incar_fo.readlines()


@pytest.fixture
def ref_incar_vasp2w90():
    with open(data_path('wannier', 'INCAR'), 'r') as reference_incar_wannier:
        yield reference_incar_wannier.readlines()


@pytest.fixture
def ref_win():
    data = Path(data_path('wannier90.win'))
    yield data.read_text()


@pytest.fixture()
def ref_retrieved():
    """Fixture: retrieved directory from an NSCF vasp run."""
    from aiida.plugins import DataFactory
    retrieved = DataFactory('folder')()
    retrieved.put_object_from_tree(path=data_path('basic_run'))
    return retrieved


@pytest.fixture()
def vasprun_parser(request):
    """Return an instance of VasprunParser for a reference vasprun.xml."""
    path, settings = path_file_and_settings('vasprun.xml', request.param)
    with open(path, 'rb') as handler:
        parser = VasprunParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def vasprun_parser_v621(request):
    """Return an instance of VasprunParser for a reference vasprun.xml of VASP6."""
    path, settings = path_file_and_settings('vasprun621.xml', request.param)
    with open(path, 'rb') as handler:
        parser = VasprunParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def neb_outcar_parser(request):
    """Return an instance of OutcarParser for a reference OUTCAR."""
    path, settings = path_file_and_settings('OUTCAR', request.param)
    with open(path, 'r') as handler:
        parser = VtstNebOutcarParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def outcar_parser(request):
    """Return an instance of OutcarParser for a reference OUTCAR."""
    path, settings = path_file_and_settings('OUTCAR', request.param)
    with open(path, 'r') as handler:
        parser = OutcarParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def poscar_parser(request):
    """Return an instance of PoscarParser for a reference POSCAR."""
    path, _ = path_file_and_settings('POSCAR', request.param)
    with open(path, 'r') as handler:
        parser = PoscarParser(handler=handler)
    return parser


@pytest.fixture()
def incar_parser(request):
    """Return an instance of IncarParser for a reference INCAR."""
    path, settings = path_file_and_settings('INCAR', request.param)
    with open(path, 'r') as handler:
        parser = IncarParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def doscar_parser(request):
    """Return an instance of DoscarParser for a reference DOSCAR."""
    path, settings = path_file_and_settings('DOSCAR', request.param)
    with open(path, 'r') as handler:
        parser = DoscarParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def chgcar_parser(request):
    """Return an instance of ChgcarParser for a reference CHGCAR."""
    path, settings = path_file_and_settings('CHGCAR', request.param)
    with open(path, 'r') as handler:
        parser = ChgcarParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def eigenval_parser(request):
    """Return an instance of EigenvalParser for a reference EIGENVAL."""
    path, settings = path_file_and_settings('EIGENVAL', request.param)
    with open(path, 'r') as handler:
        parser = EigenvalParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def kpoints_parser(request):
    """Return an instance of KpointsParser for a reference KPOINTS."""
    path, settings = path_file_and_settings('KPOINTS', request.param)
    with open(path, 'r') as handler:
        parser = KpointsParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def stream_parser(request):
    """Return an instance of StreamParser for a reference stream capture."""
    path, settings = path_file_and_settings('vasp_output', request.param)
    with open(path, 'r') as handler:
        parser = StreamParser(handler=handler, settings=settings)
    return parser


def _ref_kp_list():
    with open(data_path('kpoints', 'KPOINTS_list'), 'r') as reference_kpoints_fo:
        ref_kp_str = reference_kpoints_fo.read()
    return ref_kp_str


def _ref_kp_mesh():
    with open(data_path('kpoints', 'KPOINTS_mesh'), 'r') as reference_kpoints_fo:
        ref_kp_list = reference_kpoints_fo.read()
    return ref_kp_list


@pytest.fixture
def wannier_params():
    from aiida.orm import Dict
    return Dict(dict=dict(
        dis_num_iter=1000,
        num_bands=24,
        num_iter=0,
        num_wann=14,
        spinors=True,
    ))


@pytest.fixture
def wannier_projections():
    from aiida.orm import List
    wannier_projections = List()
    wannier_projections.extend(['Ga : s; px; py; pz', 'As : px; py; pz'])
    return wannier_projections


@pytest.fixture
def phonondb_run(tmp_path):
    phonondb = Path(data_path('phonondb'))
    copytree(phonondb, tmp_path)
    yield tmp_path


@pytest.fixture
def compare_symmetries():
    return {
        'symmetrized_cell_type': {
            'static': [
                'face centered cubic supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'face centered cubic supercell.', 'face centered cubic supercell.',
                'face centered cubic supercell.'
            ],
            'dynamic': [
                'face centered cubic supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'face centered cubic supercell.', 'face centered cubic supercell.',
                'face centered cubic supercell.'
            ]
        },
        'original_cell_type': {
            'static': [
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell'
            ],
            'dynamic': [
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell'
            ]
        },
        'num_space_group_operations': {
            'static': [48, 16, 16, 16, 16, 16, 16, 4, 4, 4, 4, 4, 4, 8, 8, 48],
            'dynamic': [48, 16, 16, 16, 16, 16, 16, 4, 4, 4, 4, 4, 4, 8, 8, 48]
        }
    }
