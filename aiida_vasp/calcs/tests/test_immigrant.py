# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import, import-outside-toplevel
"""Unit tests for VaspImmigrant calculation."""
import pytest

from aiida.engine import run
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.aiida_utils import create_authinfo


@pytest.fixture
def immigrant_with_builder(fresh_aiida_env, potcar_family, phonondb_run, localhost, mock_vasp):
    """Provide process class and inputs for importing a AiiDA-external VASP run.

    The list of objects in test_data/phonondb doesn't contain POTCAR.

    """
    from aiida_vasp.calcs.immigrant import VaspImmigrant

    create_authinfo(localhost, store=True)
    potential_family = POTCAR_FAMILY_NAME
    potential_mapping = POTCAR_MAP
    builder = VaspImmigrant.get_builder_from_folder(mock_vasp,
                                                    str(phonondb_run),
                                                    potential_family=potential_family,
                                                    potential_mapping=potential_mapping)
    # builder.potential = PotcarData.get_potcars_from_structure(builder.structure, potential_family, mapping=potential_mapping)  # pylint: disable=no-member

    # Make sure clean_workdir is not done for the immigrant (we do not want to remove the imported data)
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential'}
    for input_link in expected_inputs:
        assert builder.get(input_link, None) is not None
    return builder


@pytest.mark.usefixtures('fresh_aiida_env')
@pytest.mark.usefixtures('potcar_family')
def test_immigrant_additional(phonondb_run, localhost, mock_vasp):
    """Provide process class and inputs for importing a AiiDA-external VASP run."""
    from aiida_vasp.calcs.immigrant import VaspImmigrant
    from aiida_vasp.data.potcar import PotcarData

    create_authinfo(localhost, store=True)
    inputs = VaspImmigrant.get_inputs_from_folder(mock_vasp, str(phonondb_run), use_chgcar=True, use_wavecar=True)
    potential_family = POTCAR_FAMILY_NAME
    potential_mapping = POTCAR_MAP
    inputs.potential = PotcarData.get_potcars_from_structure(inputs.structure, potential_family, mapping=potential_mapping)
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential', 'charge_density', 'wavefunctions'}
    for input_link in expected_inputs:
        assert inputs.get(input_link, None) is not None, 'input link "{}" was not set!'.format(input_link)

    result, node = run.get_node(VaspImmigrant, **inputs)
    assert node.exit_status == 0

    # We should not have any POTCAR here
    expected_objects = ['CONTCAR', 'DOSCAR', 'EIGENVAL', 'OUTCAR', 'vasprun.xml']
    retrieved_objects = result['retrieved'].list_object_names()
    assert set(expected_objects) == set(retrieved_objects)


def test_vasp_immigrant(immigrant_with_builder):
    """Test importing a calculation from the folder of a completed VASP run."""
    builder = immigrant_with_builder

    # We need to set the parser explicitly
    # builder.metadata['options']['parser_name'] = 'vasp.vasp'
    result, node = run.get_node(builder)
    assert node.exit_status == 0

    expected_output_nodes = {'misc', 'remote_folder', 'retrieved'}
    assert expected_output_nodes.issubset(set(result))


@pytest.fixture
def immigrant_with_builder_example_3(fresh_aiida_env, potcar_family, phonondb_run, localhost, mock_vasp):  # pylint: disable=invalid-name
    """Provide process class and inputs for importing a AiiDA-external VASP run. This will be obsolete at v3."""
    from aiida_vasp.calcs.vasp import VaspCalculation

    create_authinfo(localhost, store=True)
    potential_family = POTCAR_FAMILY_NAME
    potential_mapping = POTCAR_MAP
    proc, builder = VaspCalculation.immigrant(mock_vasp,
                                              phonondb_run,
                                              potential_family=potential_family,
                                              potential_mapping=potential_mapping)
    # Make sure clean_workdir is not done for the immigrant (we do not want to remove the imported data)
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential'}
    for input_link in expected_inputs:
        assert builder.get(input_link, None) is not None
    return proc, builder


@pytest.mark.usefixtures('fresh_aiida_env')
@pytest.mark.usefixtures('potcar_family')
def test_immigrant_additional_example_3(phonondb_run, localhost, mock_vasp):  # pylint: disable=invalid-name
    """Provide process class and inputs for importing a AiiDA-external VASP run. This will be obsolete at v3."""
    from aiida_vasp.calcs.vasp import VaspCalculation
    create_authinfo(localhost, store=True)
    proc, builder = VaspCalculation.immigrant(code=mock_vasp,
                                              remote_path=phonondb_run,
                                              potential_family=POTCAR_FAMILY_NAME,
                                              potential_mapping=POTCAR_MAP,
                                              use_chgcar=True,
                                              use_wavecar=True)
    expected_inputs = {'parameters', 'structure', 'kpoints', 'potential', 'charge_density', 'wavefunctions'}
    for input_link in expected_inputs:
        assert builder.get(input_link, None) is not None, 'input link "{}" was not set!'.format(input_link)

    result, node = run.get_node(proc, **builder)
    assert node.exit_status == 0

    # We should not have any POTCAR here
    expected_objects = ['CONTCAR', 'DOSCAR', 'EIGENVAL', 'OUTCAR', 'vasprun.xml']
    retrieved_objects = result['retrieved'].list_object_names()
    assert set(expected_objects) == set(retrieved_objects)


def test_vasp_immigrant_example_3(immigrant_with_builder_example_3):  # pylint: disable=invalid-name
    """Test importing a calculation from the folder of a completed VASP run. This will be obsolete at v3."""
    immigrant, inputs = immigrant_with_builder_example_3

    # We need to set the parser explicitly
    inputs.metadata['options']['parser_name'] = 'vasp.vasp'
    result, node = run.get_node(immigrant, **inputs)
    assert node.exit_status == 0

    expected_output_nodes = {'misc', 'remote_folder', 'retrieved'}
    assert expected_output_nodes.issubset(set(result))
