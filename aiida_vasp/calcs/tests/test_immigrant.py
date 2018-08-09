# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
"""Unit tests for VaspImmigrant calculation."""
import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.aiida_utils import create_authinfo


def test_create_inputs(fresh_aiida_env, potcar_family, phonondb_run, localhost):
    """Test creating the inputs for an immigrant calculation from a folder."""
    from aiida_vasp.calcs.immigrant import VaspImmigrant
    localhost.store()
    create_authinfo(localhost, store=True)
    immigrant = VaspImmigrant(remote_workdir=phonondb_run.strpath)
    potcar_spec = {'family': POTCAR_FAMILY_NAME, 'map': POTCAR_MAP}
    with localhost.get_transport() as open_transport:
        immigrant.create_input_nodes(open_transport=open_transport, potcar_spec=potcar_spec)
    input_names = set(immigrant.get_inputs_dict())
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential_P', 'potential_S', 'potential_Zn'}
    assert expected_inputs.issubset(input_names)
