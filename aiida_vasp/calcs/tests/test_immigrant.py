# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
"""Unit tests for VaspImmigrant calculation."""
import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.aiida_utils import create_authinfo, cmp_get_transport


@pytest.fixture
def immigrant_with_builder(fresh_aiida_env, potcar_family, phonondb_run, localhost, mock_vasp):
    """Provide process class and inputs for importing a AiiDA-external VASP run."""
    from aiida_vasp.calcs.immigrant import get_immigrant_with_builder
    create_authinfo(localhost, store=True)
    potcar_spec = {'family': POTCAR_FAMILY_NAME, 'map': POTCAR_MAP}
    proc, builder = get_immigrant_with_builder(code=mock_vasp, remote_path=phonondb_run, potcar_spec=potcar_spec)
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential'}
    for input_link in expected_inputs:
        assert builder.get(input_link, None) is not None
    return proc, builder


def test_vasp_immigrant(immigrant_with_builder):
    """Test importing a calculation from the folder of a completed VASP run."""
    from aiida import work
    immigrant, inputs = immigrant_with_builder
    runner = work.runners.new_runner()

    result = runner.run(immigrant, **inputs)

    expected_output_nodes = {'output_energies', 'output_kpoints', 'output_parameters', 'output_structure', 'remote_folder', 'retrieved'}
    assert expected_output_nodes.issubset(set(result))
