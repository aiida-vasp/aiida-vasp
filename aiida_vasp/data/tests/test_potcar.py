"""Unit test the POTCAR AiiDA data structures."""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import pytest

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
