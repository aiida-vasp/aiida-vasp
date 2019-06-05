# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
"""Unit tests for VaspImmigrant calculation."""
import pytest

from aiida.engine import run
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.aiida_utils import create_authinfo, cmp_get_transport, aiida_version, cmp_version


#@pytest.mark.skip(reason="Waiting for the immigrant being moved to aiida_core")
@pytest.fixture
def immigrant_with_builder(fresh_aiida_env, potcar_family, phonondb_run, localhost, mock_vasp):
    """Provide process class and inputs for importing a AiiDA-external VASP run."""
    from aiida_vasp.calcs.vasp import VaspCalculation

    create_authinfo(localhost, store=True)
    potcar_spec = {'family': POTCAR_FAMILY_NAME, 'map': POTCAR_MAP}
    proc, builder = VaspCalculation.immigrant(code=mock_vasp, remote_path=phonondb_run, potcar_spec=potcar_spec)
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential'}
    for input_link in expected_inputs:
        assert builder.get(input_link, None) is not None
    return proc, builder


@pytest.mark.xfail(reason="Removing POTCAR from calc.raw_input not implemented yet.")
def test_immigrant_additional(fresh_aiida_env, potcar_family, phonondb_run, localhost, mock_vasp):
    """Provide process class and inputs for importing a AiiDA-external VASP run."""
    from aiida_vasp.calcs.vasp import VaspCalculation
    create_authinfo(localhost, store=True)
    potcar_spec = {'family': POTCAR_FAMILY_NAME, 'map': POTCAR_MAP}
    proc, builder = VaspCalculation.immigrant(
        code=mock_vasp, remote_path=phonondb_run, potcar_spec=potcar_spec, use_chgcar=True, use_wavecar=True)
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential', 'charge_density', 'wavefunctions'}
    for input_link in expected_inputs:
        assert builder.get(input_link, None) is not None, 'input link "{}" was not set!'.format(input_link)

    result, node = run.get_node(proc, **builder)

    assert node.exit_status == 0

    expected_files = {'INCAR', 'POSCAR', 'KPOINTS', 'CHGCAR', 'WAVECAR'}
    banned_files = {'POTCAR'}

    calc = result['retrieved']


    assert 'raw_input' in calc.folder.get_content_list()
    input_folder = calc.folder.get_subfolder('raw_input')

    input_files = set(input_folder.get_content_list())
    assert expected_files.issubset(input_files)
    assert banned_files.isdisjoint(input_files)


def test_vasp_immigrant(immigrant_with_builder):
    """Test importing a calculation from the folder of a completed VASP run."""
    immigrant, inputs = immigrant_with_builder

    # We need to set the parser explicitly
    inputs.metadata['options']['parser_name'] = 'vasp.vasp'
    result, node = run.get_node(immigrant, **inputs)
    assert node.exit_status == 0

    expected_output_nodes = {'parameters', 'remote_folder', 'retrieved'}
    assert expected_output_nodes.issubset(set(result))
