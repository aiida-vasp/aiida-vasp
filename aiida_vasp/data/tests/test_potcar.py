"""Unit test the POTCAR AiiDA data structures."""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import tarfile
import subprocess as sp
from pathlib import Path
import pytest

from pymatgen.io.vasp import PotcarSingle
from aiida.common.exceptions import UniquenessError, NotExistent

from aiida_vasp.utils.aiida_utils import get_data_node, get_data_class
from aiida_vasp.utils.fixtures.testdata import data_path, read_file
from aiida_vasp.utils.fixtures.environment import fresh_aiida_env
from aiida_vasp.utils.fixtures.data import potcar_node_pair, potcar_family, temp_pot_folder, POTCAR_MAP, legacy_potcar_family
from aiida_vasp.commands.potcar import potcar

from aiida_vasp.data.potcar import PotcarGroup, migrate_potcar_group


def test_creation(fresh_aiida_env, potcar_node_pair):
    """Test creating a data node pair."""
    potcar_node = get_data_class('vasp.potcar').find_one(symbol='As')
    file_node = potcar_node.find_file_node()
    assert potcar_node.pk == potcar_node_pair['potcar'].pk
    assert file_node.pk == potcar_node_pair['file'].pk


def test_hashing(fresh_aiida_env):
    """Ensure the file and content sha512 hash equivalently for the same POTCAR."""
    potcar_file_cls = get_data_class('vasp.potcar_file')
    potcar_path = ['potcar', 'As', 'POTCAR']

    file_sha512 = potcar_file_cls.get_file_sha512(data_path(*potcar_path))
    content_sha512 = potcar_file_cls.get_contents_sha512(read_file(*potcar_path, mode='rb'))
    assert file_sha512 == content_sha512


# pylint: disable=protected-access
def test_store_duplicate(fresh_aiida_env, potcar_node_pair):
    """
    Storing a duplicate POTCAR node must fail.

    Uniqueness constraints to test for:

        * ``sha512`` attribute must be unique
        * the combination of all other attributes must be unique
    """
    potcar_path = data_path('potcar', 'As', 'POTCAR')

    file_node = get_data_node('vasp.potcar_file', file=potcar_path)
    file_node.set_attribute('sha512', 'foo')
    with pytest.raises(UniquenessError):
        file_node.store()

    file_node = get_data_node('vasp.potcar_file', file=potcar_path)
    file_node.set_attribute('symbol', 'Ta')
    with pytest.raises(UniquenessError):
        file_node.store()

    data_node = get_data_node('vasp.potcar', potcar_file_node=potcar_node_pair['file'])
    data_node.set_attribute('sha512', 'foo')
    with pytest.raises(UniquenessError):
        data_node.store()

    data_node = get_data_node('vasp.potcar', potcar_file_node=potcar_node_pair['file'])
    data_node.set_attribute('symbol', 'Ta')
    with pytest.raises(UniquenessError):
        data_node.store()

    assert get_data_class('vasp.potcar').find_one(symbol='As')
    assert get_data_class('vasp.potcar_file').find_one(symbol='As')


def test_export_import(fresh_aiida_env, potcar_node_pair, tmp_path):
    """Exporting and importing back may not store duplicates."""
    tempfile = tmp_path / 'potcar.aiida'

    sp.call(['verdi', 'export', 'create',
             '-n', str(potcar_node_pair['file'].pk),
             '-n', str(potcar_node_pair['potcar'].pk), str(tempfile)])  # yapf: disable

    # import with same uuid
    sp.call(['verdi', 'import', str(tempfile)])
    assert get_data_class('vasp.potcar').find_one(symbol='As')
    assert get_data_class('vasp.potcar_file').find_one(symbol='As')

    # import with different uuid
    sp.call(['verdi', 'import', data_path('potcar', 'export.aiida')])
    assert get_data_class('vasp.potcar').find_one(symbol='As')
    assert get_data_class('vasp.potcar_file').find_one(symbol='As')


def test_exists(fresh_aiida_env, potcar_node_pair):
    assert get_data_class('vasp.potcar_file').exists(element='As')
    assert not get_data_class('vasp.potcar').exists(element='Xe')


def test_find(fresh_aiida_env, potcar_node_pair):
    assert get_data_class('vasp.potcar').find_one(element='As').uuid == potcar_node_pair['potcar'].uuid
    with pytest.raises(NotExistent):
        _ = get_data_class('vasp.potcar_file').find_one(element='Xe')


def test_file_get_content(fresh_aiida_env, potcar_node_pair):
    file_node_as = potcar_node_pair['file']
    original_file = Path(data_path(file_node_as.original_file_name))
    assert original_file.read_text() == file_node_as.get_content().decode()


#def test_file_get_pymatgen(fresh_aiida_env, potcar_node_pair):
#    """
#    Create a pymatgen ``PotcarSingle`` instance from a ``PotcarFileData`` node.
#
#    Test equality and completeness of the resulting object.
#    """
#    file_node_as = potcar_node_pair['file']
#    potcar_single_as = file_node_as.get_pymatgen()
#
#    assert isinstance(potcar_single_as, PotcarSingle)
#    assert file_node_as.title == potcar_single_as.keywords['TITEL']
#
#    assert potcar_single_as.data == file_node_as.get_content()


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
    assert potcar_file_cls.exists(sha512=node_file_in.sha512)


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
    assert potcar_cls.exists(sha512=node_potcar_in.sha512)


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
    indium_2.append_atom(position=[1, 0, 0], symbols='In', name='In1')
    in2_potcars = get_data_class('vasp.potcar').get_potcars_from_structure(indium_2, potcar_family, mapping={'In': 'In_d', 'In1': 'In_d'})
    assert set(in2_potcars.keys()) == {'In', 'In1'}
    in_d_potcar = get_data_class('vasp.potcar').find(family_name=potcar_family, full_name='In_d')[0]
    assert in2_potcars['In'].uuid == in_d_potcar.uuid == in2_potcars['In1'].uuid


def test_upload(fresh_aiida_env, temp_pot_folder):
    """Test uploading a family of POTCAR files."""
    family_name = 'test_family'
    family_desc = 'Test Family'
    potcar_cls = get_data_class('vasp.potcar')
    pot_dir = str(temp_pot_folder)
    potcar_ga = Path(data_path('potcar')) / 'Ga'
    assert not potcar_ga.exists()

    potcar_cls.upload_potcar_family(pot_dir, family_name, family_desc)

    assert potcar_cls.exists(element='In')
    assert potcar_cls.exists(element='As')
    assert potcar_cls.exists(element='Ga')
    assert not potcar_ga.exists()

    # this is supposed to return only one group, however it returns 8 (= number of uploaded files)
    assert [g.label for g in potcar_cls.get_potcar_groups()] == [family_name]

    assert len(potcar_cls.get_potcar_group(family_name).nodes) >= 3

    with pytest.raises(ValueError):
        potcar_cls.upload_potcar_family(pot_dir, family_name, stop_if_existing=True)
    assert not potcar_ga.exists()

    num_files, num_added, num_uploaded = potcar_cls.upload_potcar_family(pot_dir, family_name + '_new', family_desc, stop_if_existing=False)
    assert num_files >= 3
    assert num_added >= 3
    assert num_uploaded == 0
    assert not potcar_ga.exists()


def test_export_family_folder(fresh_aiida_env, potcar_family, tmp_path):
    """Test exporting to folder."""
    export_dir = tmp_path / 'export'
    export_dir.mkdir()
    potcar_cls = get_data_class('vasp.potcar')

    # Check that dry run works and does not leave anything in the directory
    potcar_cls.export_family_folder(potcar_family, path=export_dir, dry_run=True)
    elements = []
    for item in export_dir.iterdir():
        elements.append(item)
    assert not elements

    # Start check for actual export
    files = potcar_cls.export_family_folder(potcar_family, path=export_dir, dry_run=False)
    family_dir = export_dir / potcar_family
    subdirs = set(str(subpath.name) for subpath in family_dir.iterdir())
    assert set(['Ga', 'As', 'In_d']).issubset(subdirs)
    potcar_ga = family_dir / 'Ga' / 'POTCAR'
    assert potcar_ga in files
    assert len(files) >= 3
    assert potcar_ga.is_file()
    assert 'TITEL' in potcar_ga.read_text()

    new_dir = export_dir / 'new_dir'
    potcar_cls.export_family_folder(potcar_family, path=new_dir, dry_run=False)
    assert new_dir.exists()


def test_export_family_archive(fresh_aiida_env, potcar_family, tmp_path):
    """Test exporting to archive."""
    export_dir = tmp_path / 'export'
    export_dir.mkdir()
    potcar_cls = get_data_class('vasp.potcar')

    # Check that dry run works and does not leave anything in the directory
    potcar_cls.export_family_archive(potcar_family, path=export_dir, dry_run=True)
    elements = []
    for item in export_dir.iterdir():
        elements.append(item)
    assert not elements

    # Start check for actual export
    ar_path, _ = potcar_cls.export_family_archive(potcar_family, path=export_dir, dry_run=False)
    archive = tarfile.open(str(ar_path))
    assert set(['Ga/POTCAR', 'As/POTCAR', 'In_d/POTCAR']).issubset(set(archive.getnames()))
    potcar_in = archive.extractfile('In_d/POTCAR')
    try:
        content = potcar_in.read()
        assert b'TITEL' in content
    finally:
        potcar_in.close()


def test_create_equivalence(potcar_family):
    """Create from file (during upload) and from contents and ensure equivalence."""
    potcar_file_cls = get_data_class('vasp.potcar_file')
    potcar_path = ['potcar', 'As', 'POTCAR']
    potcar_file, created = potcar_file_cls.get_or_create_from_contents(read_file(*potcar_path, mode='rb'))
    assert not created
    assert potcar_file.sha512 == potcar_file_cls.find_one(element='As').sha512
    assert potcar_file.uuid == potcar_file_cls.find_one(element='As').uuid

    potcar_cls = get_data_class('vasp.potcar')
    potcar, created = potcar_cls.get_or_create_from_contents(read_file(*potcar_path, mode='rb'))
    assert not created
    assert potcar.sha512 == potcar_cls.find_one(element='As').sha512
    assert potcar.uuid == potcar_cls.find_one(element='As').uuid


def test_get_poctcars_dict(potcar_family):
    """Test the keys are the same as the input element names."""
    potcar_cls = get_data_class('vasp.potcar')
    elements = POTCAR_MAP.keys()
    mapping = POTCAR_MAP
    potcar_dict = potcar_cls.get_potcars_dict(elements=elements, family_name=potcar_family, mapping=mapping)
    assert set(potcar_dict.keys()) == set(elements)
    assert [potcar_dict[element].full_name for element in elements] == [mapping[element] for element in elements]


def test_family_migrate(potcar_family, legacy_potcar_family):
    """Test the migration from OLD potcar family to the new ones"""

    old_family_name, legacy_group_class = legacy_potcar_family
    legacy_group = legacy_group_class.objects.get(label=old_family_name)
    migrate_potcar_group()

    # Old group should still be there
    assert legacy_group_class.objects.get(label=old_family_name)

    migrated = PotcarGroup.objects.get(label=old_family_name)
    uuids_original = {node.uuid for node in legacy_group.nodes}
    uuids_migrated = {node.uuid for node in migrated.nodes}
    assert uuids_migrated == uuids_original


def test_old_style_detect(potcar_family, legacy_potcar_family):
    """Test the assestion that the potcars are found old in the legacy group not the new"""
    potcar_cls = get_data_class('vasp.potcar')
    elements = POTCAR_MAP.keys()
    mapping = POTCAR_MAP
    new_group = PotcarGroup.objects.get(label=potcar_family)
    new_group.label += '_'

    # Change the name of the legacy group to the new one so it will be matched
    legacy_group_label, legacy_group_class = legacy_potcar_family
    legacy_group = legacy_group_class.objects.get(label=legacy_group_label)
    legacy_group.label = potcar_family

    # The raise Value Error should contain hints for migrate
    with pytest.raises(NotExistent, match=r'.*found in a legacy group.*'):
        potcar_dict = potcar_cls.get_potcars_dict(elements=elements, family_name=potcar_family, mapping=mapping, auto_migrate=False)

    # Change the name back and the test should now pass
    new_group.label = new_group.label[:-1]
    potcar_dict = potcar_cls.get_potcars_dict(elements=elements, family_name=potcar_family, mapping=mapping, auto_migrate=False)
    assert set(potcar_dict.keys()) == set(elements)
    assert [potcar_dict[element].full_name for element in elements] == [mapping[element] for element in elements]

    # Test the auto-migration logic
    # Change the name again, so the only the old group matches
    new_group.label += '_'
    with pytest.raises(NotExistent):
        PotcarGroup.objects.get(label=potcar_family)

    # but as long as we do auto migrate it would be fine
    potcar_dict = potcar_cls.get_potcars_dict(elements=elements, family_name=potcar_family, mapping=mapping, auto_migrate=True)
    assert set(potcar_dict.keys()) == set(elements)
    assert [potcar_dict[element].full_name for element in elements] == [mapping[element] for element in elements]
    # Validate the migrate group
    migrated = PotcarGroup.objects.get(label=potcar_family)
    uuids_original = {node.uuid for node in legacy_group.nodes}
    uuids_migrated = {node.uuid for node in migrated.nodes}
    assert uuids_migrated == uuids_original
