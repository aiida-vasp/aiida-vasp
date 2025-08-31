"""Unit tests for importing existing VASP calculation."""

import json

import pytest
from aiida import orm
from aiida.engine import run_get_node
from click.testing import CliRunner

from aiida_vasp.calcs.immigrant import VaspCalcImporter
from aiida_vasp.commands.immigrant import import_calc


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


def test_cli_import_calc(
    aiida_profile_clean, upload_potcar, phonondb_run, localhost, mock_vasp, potcar_family_name, potcar_mapping
):
    """Test the CLI interface for importing calculations."""
    runner = CliRunner()

    # Create the potential mapping as a JSON string
    potential_mapping_json = json.dumps(potcar_mapping)

    # Test the CLI command
    result = runner.invoke(
        import_calc,
        [
            '--path',
            str(phonondb_run),
            '--code',
            mock_vasp.label,
            '--potential-family',
            potcar_family_name,
            '--potential-mapping',
            potential_mapping_json,
            '--label',
            'test-import-cli',
            '--description',
            'Test calculation imported via CLI',
            '--yes',  # Auto-confirm
            '--quiet',  # Suppress output
        ],
    )

    # Check that the command succeeded
    assert result.exit_code == 0, f'CLI command failed with output: {result.output}'

    # Check that the calculation was actually imported
    # Look for the calculation with the label we set
    qb = orm.QueryBuilder()
    qb.append(orm.CalcJobNode, filters={'label': 'test-import-cli'})
    results = qb.all()

    assert len(results) == 1, 'Expected exactly one imported calculation'
    calc_node = results[0][0]

    # Verify the calculation has the expected properties
    assert calc_node.label == 'test-import-cli'
    assert calc_node.description == 'Test calculation imported via CLI'
    assert calc_node.process_state.value == 'finished'

    # Verify the calculation has the expected inputs
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential'}
    for input_link in expected_inputs:
        assert input_link in calc_node.inputs, f'Missing input: {input_link}'

    # Verify the calculation has outputs
    assert calc_node.is_finished
    assert 'misc' in calc_node.outputs
    assert 'retrieved' in calc_node.outputs


def test_cli_import_calc_with_additional_files(
    aiida_profile_clean, upload_potcar, phonondb_run, localhost, mock_vasp, potcar_family_name, potcar_mapping
):
    """Test the CLI interface for importing calculations with additional files."""
    runner = CliRunner()

    # Create the potential mapping as a JSON string
    potential_mapping_json = json.dumps(potcar_mapping)

    # Test the CLI command with additional files
    result = runner.invoke(
        import_calc,
        [
            '--path',
            str(phonondb_run),
            '--code',
            mock_vasp.label,
            '--potential-family',
            potcar_family_name,
            '--potential-mapping',
            potential_mapping_json,
            '--include-wavecar',
            '--include-chgcar',
            '--label',
            'test-import-cli-additional',
            '--yes',  # Auto-confirm
            '--quiet',  # Suppress output
        ],
    )

    # Check that the command succeeded
    assert result.exit_code == 0, f'CLI command failed with output: {result.output}'

    # Check that the calculation was actually imported
    qb = orm.QueryBuilder()
    qb.append(orm.CalcJobNode, filters={'label': 'test-import-cli-additional'})
    results = qb.all()

    assert len(results) == 1, 'Expected exactly one imported calculation'
    calc_node = results[0][0]

    # Verify the calculation has the additional input files
    assert 'wavefunctions' in calc_node.inputs, 'Missing WAVECAR input'
    assert 'charge_density' in calc_node.inputs, 'Missing CHGCAR input'
