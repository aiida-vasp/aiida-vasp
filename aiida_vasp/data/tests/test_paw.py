"""Unit tests for the PawData AiiDA data structure"""
import pytest

from aiida_vasp.io.pymatgen_aiida.vasprun import get_data_node
from aiida_vasp.utils.fixtures.testdata import data_path


@pytest.mark.wip
def test_from_folder():
    """Create a paw data node from a directory"""
    paw_data_cls = get_data_node('vasp.paw').__class__
    paw_data = paw_data_cls.from_folder(data_path('..', 'backendtests', 'LDA', 'As'))
    paw_data.store()
    assert paw_data_cls.load_paw(element='As')[0]
