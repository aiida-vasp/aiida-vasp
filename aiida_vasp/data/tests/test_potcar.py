"""Unit test the POTCAR AiiDA data structures."""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import pytest
from py import path as py_path  # pylint: disable=no-member,no-name-in-module
from pymatgen.io.vasp import PotcarSingle
from aiida.common.exceptions import UniquenessError, NotExistent
try:
    import subprocess32 as sp
except ImportError:
    import subprocess as sp

from aiida_vasp.io.pymatgen_aiida.vasprun import get_data_node, get_data_class
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.fixtures.environment import aiida_env, fresh_aiida_env


@pytest.fixture
def potcar_node_pair(fresh_aiida_env):
    """Create a POTCAR node pair."""
    potcar_path = data_path('potcar', 'As', 'POTCAR')
    potcar_file_node = get_data_node('vasp.potcar_file', file=potcar_path)
    potcar_file_node.store()
    return {'file': potcar_file_node, 'potcar': get_data_class('vasp.potcar').find(symbol='As')}


@pytest.fixture
def potcar_family(fresh_aiida_env):
    """Create a POTCAR family."""
    family_name = 'fixture_family'
    get_data_class('vasp.potcar').upload_potcar_family(data_path('potcar'), family_name)
    return family_name


def test_creation(potcar_node_pair):
    """Test creating a data node pair."""
    potcar_node = get_data_class('vasp.potcar').find(symbol='As')
    file_node = potcar_node.find_file_node()
    assert potcar_node.pk == potcar_node_pair['potcar'].pk
    assert file_node.pk == potcar_node_pair['file'].pk


# pylint: disable=protected-access
def test_store_duplicate(potcar_node_pair):
    """
    Storing a duplicate POTCAR node must fail.

    Uniqueness constraints to test for:

        * ``md5`` attribute must be unique
        * the combination of all other attributes must be unique
    """
    potcar_path = data_path('potcar', 'As', 'POTCAR')

    file_node = get_data_node('vasp.potcar_file', file=potcar_path)
    file_node._set_attr('md5', 'foo')
    with pytest.raises(UniquenessError):
        file_node.store()

    file_node = get_data_node('vasp.potcar_file', file=potcar_path)
    file_node._set_attr('symbol', 'Ta')
    with pytest.raises(UniquenessError):
        file_node.store()

    data_node = get_data_node('vasp.potcar', potcar_file_node=potcar_node_pair['file'])
    data_node._set_attr('md5', 'foo')
    with pytest.raises(UniquenessError):
        data_node.store()

    data_node = get_data_node('vasp.potcar', potcar_file_node=potcar_node_pair['file'])
    data_node._set_attr('symbol', 'Ta')
    with pytest.raises(UniquenessError):
        data_node.store()

    assert get_data_class('vasp.potcar').find(symbol='As')
    assert get_data_class('vasp.potcar_file').find(symbol='As')


def test_export_import(potcar_node_pair, tmpdir):
    """Exporting and importing back may not store duplicates."""
    tempfile = tmpdir.join('potcar.aiida')

    sp.call(['verdi', 'export', 'create',
             '-n', str(potcar_node_pair['file'].pk),
             '-n', str(potcar_node_pair['potcar'].pk), str(tempfile)])  # yapf: disable

    # import with same uuid
    sp.call(['verdi', 'import', str(tempfile)])
    assert get_data_class('vasp.potcar').find(symbol='As')
    assert get_data_class('vasp.potcar_file').find(symbol='As')

    # import with different uuid
    sp.call(['verdi', 'import', data_path('potcar', 'export.aiida')])
    assert get_data_class('vasp.potcar').find(symbol='As')
    assert get_data_class('vasp.potcar_file').find(symbol='As')


def test_exists(potcar_node_pair):
    assert get_data_class('vasp.potcar_file').exists(element='As')
    assert not get_data_class('vasp.potcar').exists(element='Xe')


def test_find(potcar_node_pair):
    assert get_data_class('vasp.potcar').find(element='As').uuid == potcar_node_pair['potcar'].uuid
    with pytest.raises(NotExistent):
        _ = get_data_class('vasp.potcar_file').find(element='Xe')


def test_file_get_content(potcar_node_pair):
    file_node_as = potcar_node_pair['file']
    original_file = py_path.local(file_node_as.original_file_name)
    assert original_file.read() == file_node_as.get_content()


def test_file_get_pymatgen(potcar_node_pair):
    """
    Create a pymatgen ``PotcarSingle`` instance from a ``PotcarFileData`` node.

    Test equality and completeness of the resulting object.
    """
    file_node_as = potcar_node_pair['file']
    potcar_single_as = file_node_as.get_pymatgen()

    assert isinstance(potcar_single_as, PotcarSingle)
    assert file_node_as.title == potcar_single_as.keywords['TITEL']

    assert potcar_single_as.data == file_node_as.get_content()


def test_file_get_or_create(potcar_node_pair):
    """Test get_or_create of PotcarFileData."""
    potcar_as_path = data_path('potcar', 'As', 'POTCAR')
    potcar_file_cls = get_data_class('vasp.potcar_file')
    file_as = potcar_node_pair['file']
    node_file_as, created_file_as = potcar_file_cls.get_or_create(potcar_as_path)
    assert not created_file_as
    assert file_as.pk == node_file_as.pk

    potcar_in_path = data_path('potcar', 'In_d', 'POTCAR')
    node_file_in, created_file_in = potcar_file_cls.get_or_create(potcar_in_path)
    assert created_file_in
    assert potcar_file_cls.exists(md5=node_file_in.md5)


def test_potcar_get_or_create(potcar_node_pair):
    """Test get_or_create method of PotcarData."""
    potcar_cls = get_data_class('vasp.potcar')
    file_cls = get_data_class('vasp.potcar_file')
    file_as = potcar_node_pair['file']
    potcar_as = potcar_node_pair['potcar']
    node_potcar_as, created_potcar_as = potcar_cls.get_or_create(file_as)
    assert not created_potcar_as
    assert potcar_as.pk == node_potcar_as.pk

    potcar_in_path = data_path('potcar', 'In_d', 'POTCAR')
    node_potcar_in, created_potcar_in = potcar_cls.get_or_create(file_cls(file=potcar_in_path))
    assert created_potcar_in
    assert potcar_cls.exists(md5=node_potcar_in.md5)


def test_potcar_from_file(fresh_aiida_env):
    """Test creating a node pair from a file, creating the data node first."""
    potcar_cls = get_data_node('vasp.potcar')
    _, created = potcar_cls.get_or_create_from_file(data_path('potcar', 'As', 'POTCAR'))
    assert created
    _, created = potcar_cls.get_or_create_from_file(data_path('potcar', 'As', 'POTCAR'))
    assert not created


def test_potcar_from_structure(potcar_family):
    """Test getting POTCARS from a family for a structure."""
    indium_2 = get_data_node('structure')
    indium_2.append_atom(position=[0, 0, 0], symbols='In')
    indium_2.append_atom(position=[1, 0, 0], symbols='In')
    in2_potcars = get_data_class('vasp.potcar').get_potcars_dict(indium_2, potcar_family)
    assert [kind.name for kind in in2_potcars.keys()] == ['In']


def test_upload(fresh_aiida_env):
    """Test uploading a family of POTCAR files."""
    family_name = 'test_family'
    family_desc = 'Test Family'
    potcar_cls = get_data_class('vasp.potcar')

    potcar_cls.upload_potcar_family(data_path('potcar'), family_name, family_desc)

    assert potcar_cls.exists(element='In')
    assert potcar_cls.exists(element='As')
    assert potcar_cls.exists(element='Ga')

    assert [g.name for g in potcar_cls.get_potcar_groups()] == [family_name]
    assert len(potcar_cls.get_potcar_group(family_name).nodes) >= 3

    with pytest.raises(ValueError):
        potcar_cls.upload_potcar_family(data_path('potcar'), family_name, stop_if_existing=True)

    num_files, num_uploaded = potcar_cls.upload_potcar_family(data_path('potcar'), family_name, stop_if_existing=False)
    assert num_files >= 3
    assert num_uploaded == 0
