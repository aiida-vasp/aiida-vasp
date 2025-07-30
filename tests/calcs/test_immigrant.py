"""Unit tests for VaspImmigrant calculation."""

import pytest
from aiida.engine import run

from aiida_vasp.calcs.immigrant import VaspImmigrant
from aiida_vasp.calcs.vasp import VaspCalculation
from aiida_vasp.data.potcar import PotcarData
from aiida_vasp.utils.aiida_utils import create_authinfo


@pytest.fixture
def immigrant_with_builder(
    aiida_profile_clean, upload_potcar, phonondb_run, localhost, mock_vasp, potcar_family_name, potcar_mapping
):
    """Provide process class and inputs for importing a AiiDA-external VASP run.

    The list of objects in test_data/phonondb doesn't contain POTCAR.

    """
    create_authinfo(localhost, store=True)
    builder = VaspImmigrant.get_builder_from_folder(
        mock_vasp, str(phonondb_run), potential_family=potcar_family_name, potential_mapping=potcar_mapping
    )
    # Make sure clean_workdir is not done for the immigrant (we do not want to remove the imported data)
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential'}
    for input_link in expected_inputs:
        assert builder.get(input_link, None) is not None
    return builder


@pytest.mark.skip(reason='This immigrant not working with the new parser code')
def test_immigrant_additional(
    aiida_profile_clean, upload_potcar, phonondb_run, localhost, mock_vasp, potcar_family_name, potcar_mapping
):
    """Provide process class and inputs for importing a AiiDA-external VASP run."""
    create_authinfo(localhost, store=True)
    inputs = VaspImmigrant.get_inputs_from_folder(mock_vasp, str(phonondb_run), use_chgcar=True, use_wavecar=True)
    inputs.potential = PotcarData.get_potcars_from_structure(
        inputs.structure, potcar_family_name, mapping=potcar_mapping
    )
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential', 'charge_density', 'wavefunctions'}
    for input_link in expected_inputs:
        assert inputs.get(input_link, None) is not None, f'input link "{input_link}" was not set!'

    result, node = run.get_node(VaspImmigrant, **inputs)
    assert node.exit_status == 0

    # We should not have any POTCAR here
    expected_objects = ['CONTCAR', 'DOSCAR', 'EIGENVAL', 'OUTCAR', 'vasprun.xml']
    retrieved_objects = result['retrieved'].base.repository.list_object_names()
    assert set(expected_objects) == set(retrieved_objects)


@pytest.mark.skip(reason='This immigrant not working with the new parser code')
def test_vasp_immigrant(immigrant_with_builder):
    """Test importing a calculation from the folder of a completed VASP run."""
    builder = immigrant_with_builder

    # We need to set the parser explicitly
    # builder.metadata['options']['parser_name'] = 'vasp.vasp'
    result, node = run.get_node(builder)
    assert node.exit_status == 0

    expected_output_nodes = {'misc', 'retrieved'}
    assert expected_output_nodes.issubset(set(result))


@pytest.fixture
def immigrant_with_builder_example_3(
    aiida_profile_clean, upload_potcar, potcar_family_name, potcar_mapping, phonondb_run, localhost, mock_vasp
):
    """Provide process class and inputs for importing a AiiDA-external VASP run. This will be obsolete at v3."""

    create_authinfo(localhost, store=True)
    proc, builder = VaspCalculation.immigrant(
        mock_vasp, phonondb_run, potential_family=potcar_family_name, potential_mapping=potcar_mapping
    )
    # Make sure clean_workdir is not done for the immigrant (we do not want to remove the imported data)
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential'}
    for input_link in expected_inputs:
        assert builder.get(input_link, None) is not None
    return proc, builder


@pytest.mark.skip(reason='This immigrant not working with the new parser code')
def test_immigrant_additional_example_3(
    aiida_profile_clean, upload_potcar, phonondb_run, localhost, mock_vasp, potcar_family_name, potcar_mapping
):  # pylint: disable=invalid-name
    """Provide process class and inputs for importing a AiiDA-external VASP run. This will be obsolete at v3."""
    create_authinfo(localhost, store=True)
    proc, builder = VaspCalculation.immigrant(
        code=mock_vasp,
        remote_path=phonondb_run,
        potential_family=potcar_family_name,
        potential_mapping=potcar_mapping,
        use_chgcar=True,
        use_wavecar=True,
    )
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential', 'charge_density', 'wavefunctions'}
    for input_link in expected_inputs:
        assert builder.get(input_link, None) is not None, f'input link "{input_link}" was not set!'

    result, node = run.get_node(proc, **builder)
    assert node.exit_status == 0

    # We should not have any POTCAR here
    expected_objects = ['CONTCAR', 'DOSCAR', 'EIGENVAL', 'OUTCAR', 'vasprun.xml']
    retrieved_objects = result['retrieved'].base.repository.list_object_names()
    assert set(expected_objects) == set(retrieved_objects)


@pytest.mark.skip(reason='This immigrant not working with the new parser code')
def test_vasp_immigrant_example_3(immigrant_with_builder_example_3):  # pylint: disable=invalid-name
    """Test importing a calculation from the folder of a completed VASP run. This will be obsolete at v3."""
    immigrant, inputs = immigrant_with_builder_example_3

    # We need to set the parser explicitly
    inputs.metadata['options']['parser_name'] = 'vasp.vasp'
    result, node = run.get_node(immigrant, **inputs)
    assert node.exit_status == 0

    expected_output_nodes = {'misc', 'retrieved'}
    assert expected_output_nodes.issubset(set(result))
