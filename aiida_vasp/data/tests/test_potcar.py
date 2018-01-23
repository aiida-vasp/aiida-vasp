"""Unit test the POTCAR AiiDA data structures."""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import tarfile

import pytest
from py import path as py_path  # pylint: disable=no-member,no-name-in-module
from pymatgen.io.vasp import PotcarSingle
from aiida.common.exceptions import UniquenessError, NotExistent
try:
    import subprocess32 as sp
except ImportError:
    import subprocess as sp

from aiida_vasp.utils.aiida_utils import get_data_node, get_data_class
from aiida_vasp.utils.fixtures.testdata import data_path, read_file
from aiida_vasp.utils.fixtures.environment import aiida_env, fresh_aiida_env
from aiida_vasp.utils.fixtures.data import potcar_node_pair, potcar_family


def test_creation(fresh_aiida_env, potcar_node_pair):
    """Test creating a data node pair."""
    potcar_node = get_data_class('vasp.potcar').find(symbol='As')
    file_node = potcar_node.find_file_node()
    assert potcar_node.pk == potcar_node_pair['potcar'].pk
    assert file_node.pk == potcar_node_pair['file'].pk


def test_hashing(aiida_env):
    """Ensure the file and content md5 hash equivalently for the same POTCAR."""
    potcar_file_cls = get_data_class('vasp.potcar_file')
    potcar_path = ['potcar', 'As', 'POTCAR']

    file_md5 = potcar_file_cls.get_file_md5(data_path(*potcar_path))
    content_md5 = potcar_file_cls.get_contents_md5(read_file(*potcar_path))
    assert file_md5 == content_md5


# pylint: disable=protected-access
def test_store_duplicate(fresh_aiida_env, potcar_node_pair):
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


def test_export_import(fresh_aiida_env, potcar_node_pair, tmpdir):
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


def test_exists(fresh_aiida_env, potcar_node_pair):
    assert get_data_class('vasp.potcar_file').exists(element='As')
    assert not get_data_class('vasp.potcar').exists(element='Xe')


def test_find(fresh_aiida_env, potcar_node_pair):
    assert get_data_class('vasp.potcar').find(element='As').uuid == potcar_node_pair['potcar'].uuid
    with pytest.raises(NotExistent):
        _ = get_data_class('vasp.potcar_file').find(element='Xe')


def test_file_get_content(fresh_aiida_env, potcar_node_pair):
    file_node_as = potcar_node_pair['file']
    original_file = py_path.local(file_node_as.original_file_name)
    assert original_file.read() == file_node_as.get_content()


def test_file_get_pymatgen(fresh_aiida_env, potcar_node_pair):
    """
    Create a pymatgen ``PotcarSingle`` instance from a ``PotcarFileData`` node.

    Test equality and completeness of the resulting object.
    """
    file_node_as = potcar_node_pair['file']
    potcar_single_as = file_node_as.get_pymatgen()

    assert isinstance(potcar_single_as, PotcarSingle)
    assert file_node_as.title == potcar_single_as.keywords['TITEL']

    assert potcar_single_as.data == file_node_as.get_content()


def test_file_get_or_create(fresh_aiida_env, potcar_node_pair):
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


def test_potcar_get_or_create(fresh_aiida_env, potcar_node_pair):
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


def test_potcar_from_structure(fresh_aiida_env, potcar_family):
    """Test getting POTCARS from a family for a structure."""
    indium_2 = get_data_node('structure')
    indium_2.append_atom(position=[0, 0, 0], symbols='In')
    indium_2.append_atom(position=[1, 0, 0], symbols='In')
    in2_potcars = get_data_class('vasp.potcar').get_potcars_from_structure(indium_2, potcar_family)
    assert [kind[0] for kind in in2_potcars.keys()] == ['In']


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


def test_export_family_folder(fresh_aiida_env, potcar_family, tmpdir):
    """Test exporting to folder."""
    potcar_cls = get_data_class('vasp.potcar')

    potcar_cls.export_family_folder(potcar_family, path=tmpdir, dry_run=True)
    assert not tmpdir.listdir()

    files = potcar_cls.export_family_folder(potcar_family, path=tmpdir, dry_run=False)
    exportdir = tmpdir.join(potcar_family)
    subdirs = set(str(subpath.basename) for subpath in exportdir.listdir())
    assert set(['Ga', 'As', 'In_d']).issubset(subdirs)
    potcar_ga = exportdir.join('Ga', 'POTCAR')
    assert str(potcar_ga) in files
    assert len(files) >= 3
    assert potcar_ga.isfile()
    assert 'TITEL' in potcar_ga.read()

    new_dir = 'new_dir'
    potcar_cls.export_family_folder(potcar_family, path=tmpdir.join(new_dir), dry_run=False)
    assert tmpdir.join(new_dir).exists()


def test_export_family_archive(fresh_aiida_env, potcar_family, tmpdir):
    """Test exporting to archive."""
    potcar_cls = get_data_class('vasp.potcar')

    potcar_cls.export_family_archive(potcar_family, path=tmpdir, dry_run=True)
    assert not tmpdir.listdir()

    ar_path, _ = potcar_cls.export_family_archive(potcar_family, path=tmpdir, dry_run=False)
    archive = tarfile.open(str(ar_path))
    assert set(['Ga/POTCAR', 'As/POTCAR', 'In_d/POTCAR']).issubset(set(archive.getnames()))
    potcar_in = archive.extractfile('In_d/POTCAR')
    try:
        content = potcar_in.read()
        assert 'TITEL' in content
    finally:
        potcar_in.close()


def test_create_equivalence(potcar_family):
    """Create from file (during upload) and from contents and ensure equivalence."""
    potcar_file_cls = get_data_class('vasp.potcar_file')
    potcar_path = ['potcar', 'As', 'POTCAR']
    potcar_file, created = potcar_file_cls.get_or_create_from_contents(read_file(*potcar_path))
    assert not created
    assert potcar_file.md5 == potcar_file_cls.find(element='As').md5
    assert potcar_file.uuid == potcar_file_cls.find(element='As').uuid

    potcar_cls = get_data_class('vasp.potcar')
    potcar, created = potcar_cls.get_or_create_from_contents(read_file(*potcar_path))
    assert not created
    assert potcar.md5 == potcar_cls.find(element='As').md5
    assert potcar.uuid == potcar_cls.find(element='As').uuid
