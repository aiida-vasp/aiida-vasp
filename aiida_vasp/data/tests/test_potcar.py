"""Unit test the POTCAR AiiDA data structures."""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import pytest
from aiida.common.exceptions import UniquenessError
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
