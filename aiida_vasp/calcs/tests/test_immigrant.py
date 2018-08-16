# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
"""Unit tests for VaspImmigrant calculation."""
import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.aiida_utils import create_authinfo, cmp_get_transport


@pytest.fixture
def immigrant_with_builder(fresh_aiida_env, potcar_family, phonondb_run, localhost):
    """Test creating the inputs for an immigrant calculation from a folder."""
    from aiida_vasp.calcs.immigrant import get_immigrant_with_builder
    localhost.store()
    authinfo = create_authinfo(localhost, store=True)
    potcar_spec = {'family': POTCAR_FAMILY_NAME, 'map': POTCAR_MAP}
    proc, builder = get_immigrant_with_builder(code=localhost, remote_path=phonondb_run.strpath, potcar_spec=potcar_spec)
    input_names = set(immigrant.get_inputs_dict())
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential'}
    for input_link in expected_inputs:
        assert builder.input_link is not None
    return proc, builder


def test_vasp_immigrant(immigrant_with_builder):
    from aiida import work
    immigrant, inputs = immigrant_with_builder
    runner = work.runners.new_runner()

    result = runner.run(immigrant, **inputs)

    expected_output_nodes = {'output_energies', 'output_kpoints', 'output_parameters', 'output_structure', 'remote_folder', 'retrieved'}
    assert expected_output_nodes.issubset(set(result))
