"""pytest-style test fixtures"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name,too-many-function-args,
# pylint: disable=protected-access,abstract-class-instantiated,no-value-for-parameter,unexpected-keyword-arg
import os
from collections import OrderedDict
import subprocess as sp

import numpy
import pytest
from pymatgen.io.vasp import Poscar
from py import path as py_path  # pylint: disable=no-member,no-name-in-module

from aiida.orm import Computer, FolderData
from aiida.common.exceptions import NotExistent
from aiida_vasp.utils.aiida_utils import get_data_node, get_data_class
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.parsers.file_parsers.incar import IncarParser
from aiida_vasp.parsers.file_parsers.poscar import PoscarParser
from aiida_vasp.parsers.file_parsers.vasprun import VasprunParser
from aiida_vasp.parsers.file_parsers.outcar import OutcarParser

POTCAR_FAMILY_NAME = 'test_family'
POTCAR_MAP = {'In': 'In_sv', 'In_d': 'In_d', 'As': 'As', 'Ga': 'Ga', 'Si': 'Si', 'P': 'P', 'S': 'S', 'Zn': 'Zn'}


@pytest.fixture(scope='session')
def localhost_dir(tmpdir_factory):
    return tmpdir_factory.mktemp('localhost_work')


@pytest.fixture
def localhost(fresh_aiida_env, localhost_dir):
    """Fixture for a local computer called localhost. This is currently not in the AiiDA fixtures."""
    try:
        computer = Computer.objects.get(name='localhost')
    except NotExistent:
        computer = Computer(name='localhost', hostname='localhost', transport_type='local', scheduler_type='direct', workdir=localhost_dir.strpath).store()
    return computer


@pytest.fixture
def vasp_params(aiida_env):
    incar_io = IncarParser(data=get_data_class('dict')( dict={'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'sigma': 0.5, 'magmom': '30 * 2*0.'}))
    return incar_io.incar


@pytest.fixture
def potcar_node_pair(fresh_aiida_env):
    """Create a POTCAR node pair."""
    potcar_path = data_path('potcar', 'As', 'POTCAR')
    potcar_file_node = get_data_node('vasp.potcar_file', file=potcar_path)
    potcar_file_node.store()
    return {'file': potcar_file_node, 'potcar': get_data_class('vasp.potcar').find_one(symbol='As')}


@pytest.fixture
def temp_pot_folder(tmpdir):
    """A temporary copy of the potcar test data folder, to avoid extracting tar files inside the repo."""
    potcar_ga = py_path.local(data_path('potcar')).join('Ga')
    assert not potcar_ga.exists()
    pot_archive = py_path.local(data_path('potcar'))
    target = tmpdir.join('potentials')
    pot_archive.copy(target)
    return target


# pylint: disable=protected-access
def duplicate_potcar_data(potcar_node):
    """Create and store (and return) a duplicate of a given PotcarData node."""
    from aiida_vasp.data.potcar import temp_potcar
    file_node = get_data_node('vasp.potcar_file')
    with temp_potcar(potcar_node.get_content()) as potcar_file:
        file_node.add_file(potcar_file.strpath)
        file_node.set_attribute('md5', 'abcd')
        file_node.set_attribute('full_name', potcar_node.full_name)
        file_node.store()
    data_node, _ = get_data_class('vasp.potcar').get_or_create(file_node)
    return data_node


@pytest.fixture
def potcar_family(aiida_env, temp_pot_folder):
    """Create a POTCAR family."""
    potcar_ga = py_path.local(data_path('potcar')).join('Ga')
    family_name = POTCAR_FAMILY_NAME
    family_desc = 'A POTCAR family used as a test fixture. Contains only unusable POTCAR files.'
    potcar_cls = get_data_class('vasp.potcar')
    potcar_cls.upload_potcar_family(temp_pot_folder.strpath, family_name, family_desc, stop_if_existing=False)
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
def potentials(potcar_family):
    """Fixture for two incomplete POTPAW potentials"""
    potcar_cls = get_data_class('vasp.potcar')
    potentials = potcar_cls.get_potcars_dict(['In', 'In_d', 'As'], family_name=potcar_family, mapping=POTCAR_MAP)

    return potentials


@pytest.fixture(params=['cif', 'str'])
def vasp_structure(request, aiida_env):
    """Fixture: StructureData or CifData"""
    from aiida_vasp.utils.fixtures.testdata import data_path
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
    """Fixture: Well formed POSCAR contents"""
    aiida_structure = vasp_structure
    if isinstance(vasp_structure, get_data_class('cif')):
        ase_structure = vasp_structure.get_ase()
        aiida_structure = get_data_node('structure', ase=ase_structure)
    writer = PoscarParser(data=aiida_structure)
    return writer


@pytest.fixture(params=['mesh', 'list'])
def vasp_kpoints(request, aiida_env):
    """Fixture: (kpoints object, resulting KPOINTS)"""
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
def mock_vasp(aiida_env, localhost):
    """Points to a mock-up of a VASP executable."""
    from aiida.orm import Code
    from aiida.orm.querybuilder import QueryBuilder
    query_builder = QueryBuilder()
    query_builder.append(Code, tag='code')
    query_builder.add_filter('code', {'label': {'==': 'mock-vasp'}})
    query_results = query_builder.all()
    if query_results:
        code = query_results[0][0]
    else:
        os_env = os.environ.copy()
        if not localhost.pk:
            localhost.store()
        mock_vasp_path = sp.check_output(['which', 'mock-vasp'], env=os_env).strip()
        code = Code()
        code.label = 'mock-vasp'
        code.description = 'Mock VASP for tests'
        code.set_remote_computer_exec((localhost, mock_vasp_path))
        code.set_input_plugin_name('vasp.vasp')
        aiidapath = py_path.local(aiida_env.root_dir).join('.aiida')
        code.set_prepend_text('export AIIDA_PATH={}'.format(aiidapath))

    return code


@pytest.fixture()
def vasp_chgcar(aiida_env):
    """CHGCAR node and reference fixture"""
    from aiida.plugins import DataFactory
    from aiida_vasp.utils.fixtures.testdata import data_path
    chgcar_path = data_path('chgcar', 'CHGCAR')
    chgcar = DataFactory('vasp.chargedensity')(file=chgcar_path)
    with open(chgcar_path, 'r') as ref_chgcar_fo:
        ref_chgcar = ref_chgcar_fo.read()
    return chgcar, ref_chgcar


@pytest.fixture()
def vasp_wavecar(aiida_env):
    """WAVECAR node and reference fixture"""
    from aiida.plugins import DataFactory
    from aiida_vasp.utils.fixtures.testdata import data_path
    wavecar_path = data_path('wavecar', 'WAVECAR')
    wavecar = DataFactory('vasp.wavefun')(file=wavecar_path)
    with open(wavecar_path, 'r') as ref_wavecar_fo:
        ref_wavecar = ref_wavecar_fo.read()
    return wavecar, ref_wavecar


@pytest.fixture
def ref_incar():
    from aiida_vasp.utils.fixtures.testdata import data_path
    with open(data_path('test_relax_wc/inp', 'INCAR'), 'r') as reference_incar_fo:
        yield reference_incar_fo.read().strip()


@pytest.fixture
def ref_incar_vasp2w90():
    data = py_path.local(data_path('incar_set', 'INCAR.vasp2w90'))
    yield data.read()


@pytest.fixture
def ref_win():
    data = py_path.local(data_path('wannier90.win'))
    yield data.read()


@pytest.fixture()
def ref_retrieved():
    """Fixture: retrieved directory from an NSCF vasp run"""
    from aiida.plugins import DataFactory
    from aiida_vasp.utils.fixtures.testdata import data_path
    retrieved = DataFactory('folder')()
    retrieved.put_object_from_tree(path=data_path('basic_run'))
    return retrieved


@pytest.fixture(params=['vasprun'])
def vasprun_parser(request):
    """Return an instance of VasprunParser for a reference vasprun.xml."""
    from aiida_vasp.parsers.settings import ParserSettings
    file_name = 'vasprun.xml'
    path = data_path(request.param, file_name)
    parser = VasprunParser(file_path=path, settings=ParserSettings({}))
    return parser


@pytest.fixture()
def outcar_parser(request):
    """Return an instance of OutcarParser for a reference OUTCAR."""
    from aiida_vasp.parsers.settings import ParserSettings
    file_name = 'OUTCAR'
    path = data_path(request.param, file_name)
    parser = OutcarParser(file_path=path, settings=ParserSettings({}))
    return parser


def _ref_kp_list():
    from aiida_vasp.utils.fixtures.testdata import data_path
    with open(data_path('kpoints', 'KPOINTS_list'), 'r') as reference_kpoints_fo:
        ref_kp_str = reference_kpoints_fo.read()
    return ref_kp_str


def _ref_kp_mesh():
    from aiida_vasp.utils.fixtures.testdata import data_path
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
    from aiida.orm.nodes.base import List
    wannier_projections = List()
    wannier_projections.extend(['Ga : s; px; py; pz', 'As : px; py; pz'])
    return wannier_projections


@pytest.fixture
def phonondb_run(tmpdir):
    phonondb = py_path.local(data_path('phonondb'))
    phonondb.copy(tmpdir)
    yield tmpdir
