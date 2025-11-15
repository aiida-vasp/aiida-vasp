"""Unit tests for importing existing VASP calculation."""

import pytest
from aiida.engine import run_get_node

from aiida_vasp.calcs.immigrant import VaspCalcImporter


@pytest.fixture
def immigrant_with_builder(
    aiida_profile_clean, upload_potcar, phonondb_run, localhost, mock_vasp, potcar_family_name, potcar_mapping
):
    """Provide process class and inputs for importing a AiiDA-external VASP run.

    The list of objects in test_data/phonondb doesn't contain POTCAR.

    """
    builder = VaspCalcImporter.get_builder_from_folder(
        mock_vasp, str(phonondb_run), potential_family=potcar_family_name, potential_mapping=potcar_mapping
    )
    # Make sure clean_workdir is not done for the immigrant (we do not want to remove the imported data)
    return builder


def test_get_builder(immigrant_with_builder):
    """Test getting the builder from an existing calculation."""
    builder = immigrant_with_builder
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential'}
    for input_link in expected_inputs:
        assert builder.get(input_link, None) is not None


def test_vasp_immigrant(immigrant_with_builder):
    """Test importing a calculation from an existing folder of a completed VASP run."""
    builder = immigrant_with_builder

    # We need to set the parser explicitly
    # builder.metadata['options']['parser_name'] = 'vasp.vasp'
    result, node = run_get_node(builder)
    assert node.exit_status == 0

    expected_output_nodes = {'misc', 'retrieved'}
    assert expected_output_nodes.issubset(set(result))


def test_immigrant_additional(mock_vasp, phonondb_run, potcar_family_name, potcar_mapping):
    """Test importing additional files from a completed VASP run."""
    builder = VaspCalcImporter.get_builder_from_folder(
        mock_vasp,
        str(phonondb_run),
        include_chgcar=True,
        include_wavecar=True,
        potential_family=potcar_family_name,
        potential_mapping=potcar_mapping,
    )

    builder.settings = {'ADDITIONAL_RETRIEVE_LIST': ['DOSCAR', 'EIGENVAL']}
    result, node = run_get_node(builder)
    assert node.exit_status == 0

    # We should not have any POTCAR here
    expected_objects = ['CONTCAR', 'OUTCAR', 'vasprun.xml', 'vasp_output', 'DOSCAR', 'EIGENVAL']
    retrieved_objects = result['retrieved'].base.repository.list_object_names()
    assert set(expected_objects) == set(retrieved_objects)
