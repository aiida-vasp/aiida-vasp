"""Test the POTCAR parser."""

import pytest

from aiida_vasp.data.potcar import PotcarData, PotcarFileData
from aiida_vasp.parsers.content_parsers.potcar import MultiPotcarIo, PotcarIo


def verify_potcario(potcario):
    assert potcario.node
    assert potcario.file_node
    assert potcario.content


def test_potcar_from_path(aiida_profile, data_path):
    """Create a PotcarIo instance from a file path."""
    potcar_path_as = data_path('potcar', 'As', 'POTCAR')
    from_ctor = PotcarIo(path=potcar_path_as)
    from_from = PotcarIo.from_(potcar_path_as)
    verify_potcario(from_from)
    assert from_ctor == from_from


def test_potcar_from_file_node(upload_potcar):
    """Create a PotcarIo instance from a PotcarFileData node."""
    potcar_file_in = PotcarFileData.find_one(element='In')
    from_ctor = PotcarIo(potcar_file_node=potcar_file_in)
    verify_potcario(from_ctor)
    from_from = PotcarIo.from_(potcar_file_in)
    assert from_ctor == from_from


def test_potcar_from_node(upload_potcar):
    """Create a PotcarIo instance from a PotcarData node."""
    potcar_ga = PotcarData.find_one(element='Ga')
    from_ctor = PotcarIo(potcar_node=potcar_ga)
    verify_potcario(from_ctor)
    from_from = PotcarIo.from_(potcar_ga)
    assert from_ctor == from_from


def test_potcar_from_contents(upload_potcar, read_file):
    """Create a PotcarIo from contents of a POTCAR file."""
    contents_as = read_file('potcar', 'As', 'POTCAR')
    from_ctor = PotcarIo(contents=contents_as.encode('utf-8'))
    verify_potcario(from_ctor)
    assert from_ctor.node.uuid == PotcarData.find_one(element='As').uuid
    from_from = PotcarIo.from_(contents_as)
    assert from_ctor == from_from


def test_file_contents_equivalence(aiida_profile, data_path, read_file):
    potcar_path_as = ['potcar', 'As', 'POTCAR']
    from_file = PotcarIo(path=data_path(*potcar_path_as))
    from_contents = PotcarIo(contents=read_file(*potcar_path_as).encode('utf-8'))
    assert from_file.sha512 == from_contents.sha512


def test_multi_round_trip(upload_potcar, potcar_family_name, tmp_path, potcar_mapping):
    """Write multiple POTCAR potentials to a file and recover the nodes stored in the db."""
    test_dir = tmp_path / 'round_trip'
    test_dir.mkdir()
    potcar_cls = PotcarData
    multi = MultiPotcarIo(
        potcar_cls.get_potcars_dict(
            elements=potcar_mapping.keys(), family_name=potcar_family_name, mapping=potcar_mapping
        ).values()
    )
    tempfile = test_dir / 'POTCAR'
    multi.write(tempfile)
    recovered = multi.read(tempfile)
    uuids_start = [potcar.node.uuid for potcar in multi.potcars]
    uuids_recov = [potcar.node.uuid for potcar in recovered.potcars]
    assert uuids_start == uuids_recov


@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_multi_from_structure(upload_potcar, potcar_family_name, vasp_structure_poscar, potcar_mapping):
    potcar_cls = PotcarData
    potcar_dict = potcar_cls.get_potcars_dict(
        elements=['As', 'In', 'In_d'], family_name=potcar_family_name, mapping=potcar_mapping
    )
    multi = MultiPotcarIo.from_structure(structure=vasp_structure_poscar._content_data, potentials_dict=potcar_dict)  # pylint: disable=protected-access
    assert [potcar.node.full_name for potcar in multi.potcars] == ['In_sv', 'As', 'In_d', 'As']
