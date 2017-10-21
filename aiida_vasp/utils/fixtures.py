"""pytest-style test fixtures"""
# pylint: disable=unused-argument,redefined-outer-name
import tempfile
import shutil
from collections import OrderedDict

import numpy
import pytest

from aiida.utils.fixtures import fixture_manager


@pytest.fixture(scope='session')
def aiida_env():
    with fixture_manager() as manager:
        yield manager


@pytest.fixture(scope='function')
def fresh_aiida_env(aiida_env):
    yield
    aiida_env.reset_db()


@pytest.fixture(scope='function')
def localhost(aiida_env):
    """Fixture for a local computer called localhost"""
    from aiida.orm import Computer
    workdir = tempfile.mkdtemp()
    computer = Computer(
        name='localhost',
        description='description',
        hostname='localhost',
        workdir=workdir,
        transport_type='local',
        scheduler_type='direct',
        enabled_state=True)
    yield computer
    shutil.rmtree(workdir)


@pytest.fixture()
def vasp_params(aiida_env):
    from aiida.orm import DataFactory

    return DataFactory('parameter')(dict=OrderedDict([('gga', 'PE'), (
        'gga_compat', False), ('lorbit', 11), ('sigma', 0.5)]))


@pytest.fixture()
def paws(aiida_env):
    """Fixture for two incomplete POTPAW potentials"""
    from aiida.orm import DataFactory
    from aiida_vasp.backendtests.common import subpath
    DataFactory('vasp.paw').import_family(
        subpath('../backendtests/LDA'),
        familyname='TEST',
        family_desc='test data',
    )
    paw_nodes = {
        'In': DataFactory('vasp.paw').load_paw(element='In')[0],
        'As': DataFactory('vasp.paw').load_paw(element='As')[0]
    }
    return paw_nodes


@pytest.fixture()
def cif_structure(aiida_env):
    from aiida_vasp.backendtests.common import subpath
    from aiida.orm import DataFactory
    cif_path = subpath('data/EntryWithCollCode43360.cif')
    return DataFactory('cif').get_or_create(cif_path)[0]


@pytest.fixture()
def aiida_structure(aiida_env):
    """Fixture for an aiida StructureData"""
    from aiida.orm import DataFactory
    larray = numpy.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    alat = 6.058
    structure = DataFactory('structure')(cell=larray * alat)
    structure.append_atom(position=[0, 0, 0], symbols='In')
    structure.append_atom(position=[.25, .25, .25], symbols='As')
    return structure


@pytest.fixture()
def kpoints_mesh(aiida_env):
    from aiida.orm import DataFactory
    kpoints = DataFactory('array.kpoints')()
    kpoints.set_kpoints_mesh([2, 2, 2])
    return kpoints


@pytest.fixture()
def kpoints_list(aiida_env):
    from aiida.orm import DataFactory
    kpoints = DataFactory('array.kpoints')()
    kpoints.set_kpoints([[0., 0., 0.], [0., 0., .5]], weights=[1., 1.])
    return kpoints


@pytest.fixture()
def vasp_code(localhost, fresh_aiida_env):
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
def ref_incar():
    from aiida_vasp.backendtests.common import subpath
    with open(subpath('data/INCAR'), 'r') as reference_incar_fo:
        yield reference_incar_fo.read()
