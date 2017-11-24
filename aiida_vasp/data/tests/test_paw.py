"""Unit tests for the PawData AiiDA data structure"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name

from aiida_vasp.io.pymatgen_aiida.vasprun import get_data_node
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.fixtures.environment import aiida_env


def test_from_folder(aiida_env):
    """Create a paw data node from a directory"""
    paw_data_cls = get_data_node('vasp.paw').__class__
    paw_data = paw_data_cls.from_folder(data_path('..', 'backendtests', 'LDA', 'As'))
    paw_data.store()
    assert paw_data_cls.load_paw(element='As')[0]
