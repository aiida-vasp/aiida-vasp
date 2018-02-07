"""pytest-style test fixtures"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import os
from collections import OrderedDict

import numpy
import pytest
from pymatgen.io.vasp import Poscar
from py import path as py_path  # pylint: disable=no-member,no-name-in-module

from aiida_vasp.utils.aiida_utils import get_data_node, get_data_class
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.incar import IncarIo

POTCAR_FAMILY_NAME = 'test_family'
POTCAR_MAP = {'In': 'In_d', 'As': 'As', 'Ge': 'Ge'}


@pytest.fixture(scope='session')
def localhost_dir(tmpdir_factory):
    return tmpdir_factory.mktemp('localhost_work')


@pytest.fixture
def localhost(aiida_env, localhost_dir):
    """Fixture for a local computer called localhost"""
    from aiida.orm import Computer
    from aiida.orm.querybuilder import QueryBuilder
    query_builder = QueryBuilder()
    query_builder.append(Computer, tag='comp')
    query_builder.add_filter('comp', {'name': {'==': 'localhost'}})
    query_results = query_builder.all()
    if query_results:
        computer = query_results[0][0]
    else:
        computer = Computer(
            name='localhost',
            description='description',
            hostname='localhost',
            workdir=localhost_dir.strpath,
            transport_type='local',
            scheduler_type='direct',
            enabled_state=True)
    return computer


@pytest.fixture
def vasp_params(aiida_env):
    incar_io = IncarIo(incar_dict={'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'sigma': 0.5, 'magmom': '30 * 2*0.'})
    return incar_io.get_param_node()


@pytest.fixture
def potcar_node_pair(aiida_env):
    """Create a POTCAR node pair."""
    potcar_path = data_path('potcar', 'As', 'POTCAR')
    potcar_file_node = get_data_node('vasp.potcar_file', file=potcar_path)
    potcar_file_node.store()
    return {'file': potcar_file_node, 'potcar': get_data_class('vasp.potcar').find(symbol='As')}


@pytest.fixture
def temp_pot_folder(tmpdir):
    """A temporary copy of the potcar test data folder, to avoid extracting tar files inside the repo."""
    pot_archive = py_path.local(data_path('potcar'))
    target = tmpdir.join('potcar')
    pot_archive.copy(target)
    return target


@pytest.fixture
def potcar_family(aiida_env, temp_pot_folder):
    """Create a POTCAR family."""
    family_name = POTCAR_FAMILY_NAME
    family_desc = 'A POTCAR family used as a test fixture. Contains only unusable POTCAR files.'
    potcar_cls = get_data_class('vasp.potcar')
    potcar_cls.upload_potcar_family(temp_pot_folder.strpath, family_name, family_desc, stop_if_existing=False)
    assert 'As' in potcar_cls.get_full_names(POTCAR_FAMILY_NAME, 'As')
    assert 'Ga' in potcar_cls.get_full_names(POTCAR_FAMILY_NAME, 'Ga')
    assert 'In_d' in potcar_cls.get_full_names(POTCAR_FAMILY_NAME, 'In')
    return family_name


@pytest.fixture
def potentials(potcar_family):
    """Fixture for two incomplete POTPAW potentials"""
    potcar_cls = get_data_class('vasp.potcar')
    potentials = potcar_cls.get_potcars_dict(['In', 'As'], family_name=potcar_family, mapping=POTCAR_MAP)

    return potentials


@pytest.fixture(params=['cif', 'str'])
def vasp_structure(request, aiida_env):
    """Fixture: StructureData or CifData"""
    from aiida_vasp.backendtests.common import subpath
    from aiida.orm import DataFactory
    if request.param == 'cif':
        cif_path = subpath('data', 'EntryWithCollCode43360.cif')
        structure = DataFactory('cif').get_or_create(cif_path)[0]
    elif request.param == 'str':
        larray = numpy.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
        alat = 6.058
        structure = DataFactory('structure')(cell=larray * alat)
        structure.append_atom(position=[0, 0, 0], symbols='In')
        structure.append_atom(position=[.25, .25, .25], symbols='As')
        structure.append_atom(position=[.5, .5, .5], symbols='In')
        structure.append_atom(position=[.75, .75, .75], symbols='As')
    return structure


@pytest.fixture()
def vasp_structure_poscar(vasp_structure):
    """Fixture: Well formed POSCAR contents"""
    ase_structure = vasp_structure.get_ase()
    aiida_structure = get_data_node('structure', ase=ase_structure)
    pmg_structure = aiida_structure.get_pymatgen()
    pmg_structure.sort()
    pmg_poscar = Poscar(pmg_structure)
    return pmg_poscar


@pytest.fixture(params=['mesh', 'list'])
def vasp_kpoints(request, aiida_env):
    """Fixture: (kpoints object, resulting KPOINTS)"""
    from aiida.orm import DataFactory
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
    """Fixture for a vasp code, the executable it points to does not exist"""
    from aiida.orm import Code
    localhost.store()
    code = Code()
    code.label = 'vasp'
    code.description = 'VASP code'
    code.set_remote_computer_exec((localhost, '/usr/local/bin/vasp'))
    code.set_input_plugin_name('vasp.vasp')
    return code


@pytest.fixture()
def vasp_chgcar(aiida_env):
    """CHGCAR node and reference fixture"""
    from aiida.orm import DataFactory
    from aiida_vasp.backendtests.common import subpath
    chgcar_path = subpath('data', 'CHGCAR')
    chgcar = DataFactory('vasp.chargedensity')(file=chgcar_path)
    with open(chgcar_path, 'r') as ref_chgcar_fo:
        ref_chgcar = ref_chgcar_fo.read()
    return chgcar, ref_chgcar


@pytest.fixture()
def vasp_wavecar(aiida_env):
    """WAVECAR node and reference fixture"""
    from aiida.orm import DataFactory
    from aiida_vasp.backendtests.common import subpath
    wavecar_path = subpath('data', 'WAVECAR')
    wavecar = DataFactory('vasp.wavefun')(file=wavecar_path)
    with open(wavecar_path, 'r') as ref_wavecar_fo:
        ref_wavecar = ref_wavecar_fo.read()
    return wavecar, ref_wavecar


@pytest.fixture()
def ref_incar():
    from aiida_vasp.backendtests.common import subpath
    with open(subpath('data', 'INCAR'), 'r') as reference_incar_fo:
        yield reference_incar_fo.read().strip()


@pytest.fixture()
def ref_retrieved_nscf():
    """Fixture: retrieved directory from an NSCF vasp run"""
    from aiida.orm import DataFactory
    from aiida_vasp.backendtests.common import subpath
    retrieved = DataFactory('folder')()
    for fname in os.listdir(subpath('data', 'retrieved_nscf', 'path')):
        retrieved.add_path(subpath('data', 'retrieved_nscf', 'path', fname), '')
    return retrieved


def _ref_kp_list():
    from aiida_vasp.backendtests.common import subpath
    with open(subpath('data', 'KPOINTS.list'), 'r') as reference_kpoints_fo:
        ref_kp_str = reference_kpoints_fo.read()
    return ref_kp_str


def _ref_kp_mesh():
    from aiida_vasp.backendtests.common import subpath
    with open(subpath('data', 'KPOINTS.mesh'), 'r') as reference_kpoints_fo:
        ref_kp_list = reference_kpoints_fo.read()
    return ref_kp_list
