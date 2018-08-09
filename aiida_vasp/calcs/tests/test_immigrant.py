# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
"""Unit tests for VaspImmigrant calculation."""
import pytest

from aiida_vasp.utils.fixtures import *


def test_create_inputs(pnonondb_run, localhost):
    """Test creating the inputs for an immigrant calculation from a folder."""
    from aiida_vasp.calcs.immigrant import VaspImmigrant
    immigrant = VaspImmigrant(remote_workdir=phonondb_run.strpath)
    with localhost.get_transport() as open_transport:
        immigrant.create_input_nodes(open_transport=open_transport)
    input_names = set(immigrant.get_inputs_dict())
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential_P', 'potential_S', 'potential_Zn'}
    assert expected_inputs.issubset(input_names)
