"""Unit test the POTCAR AiiDA data structures."""

import pathlib
import subprocess as sp
import tarfile

import pytest
from aiida import orm
from aiida.common.exceptions import NotExistent, UniquenessError

from aiida_vasp.data.potcar import PotcarData, PotcarFileData, PotcarGroup, migrate_potcar_group


def test_creation(aiida_profile_clean, potcar_node_pair):
    """Test creating a data node pair."""
    potcar_node = PotcarData.find_one(symbol='As')
    file_node = potcar_node.find_file_node()
    assert potcar_node.pk == potcar_node_pair['potcar'].pk
    assert file_node.pk == potcar_node_pair['file'].pk


def test_hashing(aiida_profile_clean, data_path, read_file):
    """Ensure the file and content sha512 hash equivalently for the same POTCAR."""
    potcar_file_cls = PotcarFileData
    potcar_path = ['potcar', 'As', 'POTCAR']

    file_sha512 = potcar_file_cls.get_file_sha512(data_path(*potcar_path))
    content_sha512 = potcar_file_cls.get_contents_sha512(read_file(*potcar_path, mode='rb'))
    assert file_sha512 == content_sha512


# pylint: disable=protected-access
def test_store_duplicate(aiida_profile_clean, potcar_node_pair, data_path):
    """
    Storing a duplicate POTCAR node must fail.

    Uniqueness constraints to test for:

        * ``sha512`` attribute must be unique
        * the combination of all other attributes must be unique
    """
    potcar_path = data_path('potcar', 'As', 'POTCAR')

    file_node = PotcarFileData(file=potcar_path)
    file_node.base.attributes.set('sha512', 'foo')
    with pytest.raises(UniquenessError):
        file_node.store()

    file_node = PotcarFileData(file=potcar_path)
    file_node.base.attributes.set('symbol', 'Ta')
    with pytest.raises(UniquenessError):
        file_node.store()

    data_node = PotcarData(potcar_file_node=potcar_node_pair['file'])
    data_node.base.attributes.set('sha512', 'foo')
    with pytest.raises(UniquenessError):
        data_node.store()

    data_node = PotcarData(potcar_file_node=potcar_node_pair['file'])
    data_node.base.attributes.set('symbol', 'Ta')
    with pytest.raises(UniquenessError):
        data_node.store()

    assert PotcarData.find_one(symbol='As')
    assert PotcarFileData.find_one(symbol='As')


def test_export_import(aiida_profile_clean, potcar_node_pair, tmp_path, data_path):
    """Exporting and importing back may not store duplicates."""
    tempfile = tmp_path / 'potcar.aiida'

    sp.call(['verdi', 'export', 'create',
             '-n', str(potcar_node_pair['file'].pk),
             '-n', str(potcar_node_pair['potcar'].pk), str(tempfile)])  # yapf: disable

    # import with same uuid
    sp.call(['verdi', 'import', str(tempfile)])
    assert PotcarData.find_one(symbol='As')
    assert PotcarFileData.find_one(symbol='As')

    # import with different uuid
    sp.call(['verdi', 'import', data_path('potcar', 'export.aiida')])
    assert PotcarData.find_one(symbol='As')
    assert PotcarFileData.find_one(symbol='As')


def test_exists(aiida_profile_clean, potcar_node_pair):
    assert PotcarFileData.exists(element='As')
    assert not PotcarData.exists(element='Xe')


def test_find(aiida_profile_clean, potcar_node_pair):
    assert PotcarData.find_one(element='As').uuid == potcar_node_pair['potcar'].uuid
    with pytest.raises(NotExistent):
        _ = PotcarFileData.find_one(element='Xe')


def test_file_get_content(aiida_profile_clean, potcar_node_pair, data_path):
    file_node_as = potcar_node_pair['file']
    original_file = pathlib.Path(data_path(file_node_as.original_file_name))
    assert original_file.read_text(encoding='utf8') == file_node_as.get_content().decode()


def test_file_get_or_create(aiida_profile_clean, potcar_node_pair, data_path):
    """Test get_or_create of PotcarFileData."""
    potcar_as_path = data_path('potcar', 'As', 'POTCAR')
    potcar_file_cls = PotcarFileData
    file_as = potcar_node_pair['file']
    node_file_as, created_file_as = potcar_file_cls.get_or_create(potcar_as_path)
    assert not created_file_as
    assert file_as.pk == node_file_as.pk

    potcar_in_path = data_path('potcar', 'In_d', 'POTCAR')
    node_file_in, created_file_in = potcar_file_cls.get_or_create(potcar_in_path)
    assert created_file_in
    assert potcar_file_cls.exists(sha512=node_file_in.sha512)


def test_potcar_get_or_create(aiida_profile_clean, potcar_node_pair, data_path):
    """Test get_or_create method of PotcarData."""
    potcar_cls = PotcarData
    file_cls = PotcarFileData
    file_as = potcar_node_pair['file']
    potcar_as = potcar_node_pair['potcar']
    node_potcar_as, created_potcar_as = potcar_cls.get_or_create(file_as)
    assert not created_potcar_as
    assert potcar_as.pk == node_potcar_as.pk

    potcar_in_path = data_path('potcar', 'In_d', 'POTCAR')
    node_potcar_in, created_potcar_in = potcar_cls.get_or_create(file_cls(file=potcar_in_path))
    assert created_potcar_in
    assert potcar_cls.exists(sha512=node_potcar_in.sha512)


def test_potcar_from_file(aiida_profile_clean, data_path):
    """Test creating a node pair from a file, creating the data node first."""
    potcar_cls = PotcarData
    _, created = potcar_cls.get_or_create_from_file(data_path('potcar', 'As', 'POTCAR'))
    assert created
    _, created = potcar_cls.get_or_create_from_file(data_path('potcar', 'As', 'POTCAR'))
    assert not created


def test_potcar_from_structure(aiida_profile_clean, upload_potcar, potcar_family_name):
    """Test getting POTCARS from a family for a structure."""
    indium_2 = orm.StructureData()
    indium_2.append_atom(position=[0, 0, 0], symbols='In')
    indium_2.append_atom(position=[1, 0, 0], symbols='In', name='In1')
    in2_potcars = PotcarData.get_potcars_from_structure(
        indium_2, potcar_family_name, mapping={'In': 'In_d', 'In1': 'In_d'}
    )
    assert set(in2_potcars.keys()) == {'In', 'In1'}
    in_d_potcar = PotcarData.find(family_name=potcar_family_name, full_name='In_d')[0]
    assert in2_potcars['In'].uuid == in_d_potcar.uuid == in2_potcars['In1'].uuid


def test_upload(aiida_profile_clean, temp_pot_folder, data_path):
    """Test uploading a family of POTCAR files."""
    family_name = 'test_family'
    family_desc = 'Test Family'
    potcar_cls = PotcarData
    pot_dir = str(temp_pot_folder)
    potcar_ga = pathlib.Path(data_path('potcar')) / 'Ga'
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

    num_files, num_added, num_uploaded = potcar_cls.upload_potcar_family(
        pot_dir, family_name + '_new', family_desc, stop_if_existing=False
    )
    assert num_files >= 3
    assert num_added >= 3
    assert num_uploaded == 0
    assert not potcar_ga.exists()


def test_export_family_folder(aiida_profile_clean, upload_potcar, potcar_family_name, tmp_path):
    """Test exporting to folder."""
    export_dir = tmp_path / 'export'
    export_dir.mkdir()
    potcar_cls = PotcarData

    # Check that dry run works and does not leave anything in the directory
    potcar_cls.export_family_folder(potcar_family_name, path=export_dir, dry_run=True)
    elements = []
    for item in export_dir.iterdir():
        elements.append(item)
    assert not elements

    # Start check for actual export
    files = potcar_cls.export_family_folder(potcar_family_name, path=export_dir, dry_run=False)
    family_dir = export_dir / potcar_family_name
    subdirs = set(str(subpath.name) for subpath in family_dir.iterdir())
    assert set(['Ga', 'As', 'In_d']).issubset(subdirs)
    potcar_ga = family_dir / 'Ga' / 'POTCAR'
    assert potcar_ga in files
    assert len(files) >= 3
    assert potcar_ga.is_file()
    assert 'TITEL' in potcar_ga.read_text()

    new_dir = export_dir / 'new_dir'
    potcar_cls.export_family_folder(potcar_family_name, path=new_dir, dry_run=False)
    assert new_dir.exists()


def test_export_family_archive(aiida_profile_clean, upload_potcar, potcar_family_name, tmp_path):
    """Test exporting to archive."""
    export_dir = tmp_path / 'export'
    export_dir.mkdir()
    potcar_cls = PotcarData

    # Check that dry run works and does not leave anything in the directory
    potcar_cls.export_family_archive(potcar_family_name, path=export_dir, dry_run=True)
    elements = []
    for item in export_dir.iterdir():
        elements.append(item)
    assert not elements

    # Start check for actual export
    ar_path, _ = potcar_cls.export_family_archive(potcar_family_name, path=export_dir, dry_run=False)
    with tarfile.open(str(ar_path)) as archive:
        assert set(['Ga/POTCAR', 'As/POTCAR', 'In_d/POTCAR']).issubset(set(archive.getnames()))
        potcar_in = archive.extractfile('In_d/POTCAR')
        content = potcar_in.read()
        assert b'TITEL' in content


def test_create_equivalence(upload_potcar, read_file):
    """Create from file (during upload) and from contents and ensure equivalence."""
    potcar_file_cls = PotcarFileData
    potcar_path = ['potcar', 'As', 'POTCAR']
    potcar_file, created = potcar_file_cls.get_or_create_from_contents(read_file(*potcar_path, mode='rb'))
    assert not created
    assert potcar_file.sha512 == potcar_file_cls.find_one(element='As').sha512
    assert potcar_file.uuid == potcar_file_cls.find_one(element='As').uuid

    potcar_cls = PotcarData
    potcar, created = potcar_cls.get_or_create_from_contents(read_file(*potcar_path, mode='rb'))
    assert not created
    assert potcar.sha512 == potcar_cls.find_one(element='As').sha512
    assert potcar.uuid == potcar_cls.find_one(element='As').uuid


def test_get_poctcars_dict(upload_potcar, potcar_family_name, potcar_mapping):
    """Test the keys are the same as the input element names."""
    potcar_cls = PotcarData
    elements = potcar_mapping.keys()
    potcar_dict = potcar_cls.get_potcars_dict(elements=elements, family_name=potcar_family_name, mapping=potcar_mapping)
    assert set(potcar_dict.keys()) == set(elements)
    assert [potcar_dict[element].full_name for element in elements] == [potcar_mapping[element] for element in elements]


def test_family_migrate(upload_potcar, legacy_potcar_family):
    """Test the migration from OLD potcar family to the new ones"""

    old_family_name, legacy_group_class = legacy_potcar_family
    legacy_group = legacy_group_class.collection.get(label=old_family_name)
    migrate_potcar_group()

    # Old group should still be there
    assert legacy_group_class.collection.get(label=old_family_name)

    migrated = PotcarGroup.collection.get(label=old_family_name)
    uuids_original = {node.uuid for node in legacy_group.nodes}
    uuids_migrated = {node.uuid for node in migrated.nodes}
    assert uuids_migrated == uuids_original


def test_old_style_detect(aiida_profile_clean, potcar_family_name, potcar_mapping, legacy_potcar_family):
    """Test the assestion that the potcars are found old in the legacy group not the new"""
    potcar_cls = PotcarData
    elements = potcar_mapping.keys()
    new_group = PotcarGroup.collection.get(label=potcar_family_name)
    new_group.label += '_'

    # Change the name of the legacy group to the new one so it will be matched
    legacy_group_label, legacy_group_class = legacy_potcar_family
    legacy_group = legacy_group_class.collection.get(label=legacy_group_label)
    legacy_group.label = potcar_family_name

    # The raise Value Error should contain hints for migrate
    with pytest.raises(NotExistent, match=r'.*found in a legacy group.*'):
        potcar_dict = potcar_cls.get_potcars_dict(
            elements=elements,
            family_name=potcar_family_name,
            mapping=potcar_mapping,
            auto_migrate=False,
        )

    # Change the name back and the test should now pass
    new_group.label = new_group.label[:-1]
    potcar_dict = potcar_cls.get_potcars_dict(
        elements=elements,
        family_name=potcar_family_name,
        mapping=potcar_mapping,
        auto_migrate=False,
    )
    assert set(potcar_dict.keys()) == set(elements)
    assert [potcar_dict[element].full_name for element in elements] == [potcar_mapping[element] for element in elements]

    # Test the auto-migration logic
    # Change the name again, so the only the old group matches
    new_group.label += '_'
    with pytest.raises(NotExistent):
        PotcarGroup.collection.get(label=potcar_family_name)

    # but as long as we do auto migrate it would be fine
    potcar_dict = potcar_cls.get_potcars_dict(
        elements=elements, family_name=potcar_family_name, mapping=potcar_mapping, auto_migrate=True
    )
    assert set(potcar_dict.keys()) == set(elements)
    assert [potcar_dict[element].full_name for element in elements] == [potcar_mapping[element] for element in elements]
    # Validate the migrate group
    migrated = PotcarGroup.collection.get(label=potcar_family_name)
    uuids_original = {node.uuid for node in legacy_group.nodes}
    uuids_migrated = {node.uuid for node in migrated.nodes}
    assert uuids_migrated == uuids_original
