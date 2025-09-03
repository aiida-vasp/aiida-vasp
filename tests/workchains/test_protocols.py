import pytest

from aiida_vasp.workchains.v2 import (
    VaspBandsWorkChain,
    VaspConvergenceWorkChain,
    VaspHybridBandsWorkChain,
    VaspRelaxWorkChain,
    VaspWorkChain,
)


@pytest.fixture
def basic_env(aiida_profile, mock_vasp, potcar_family_name, upload_potcar):
    """
    Test defining the inputs for a VASP workchain
    """
    pass


@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_vasp_protocol(basic_env, mock_vasp, vasp_structure, potcar_family_name):
    """Test VASP workchain protocol"""

    builder = VaspWorkChain.get_builder_from_protocol(
        code=mock_vasp,
        structure=vasp_structure,
        overrides={'potential_family': potcar_family_name, 'potential_mapping': {'In_d': 'In_d'}},
    )

    assert builder.structure == vasp_structure
    assert builder.code == mock_vasp
    assert builder.parameters is not None


@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_relax_protocol(basic_env, mock_vasp, vasp_structure, potcar_family_name):
    """Test VASP relax workchain protocol"""

    builder = VaspRelaxWorkChain.get_builder_from_protocol(
        code=mock_vasp,
        structure=vasp_structure,
        overrides={'vasp': {'potential_family': potcar_family_name, 'potential_mapping': {'In_d': 'In_d'}}},
    )

    assert builder.structure == vasp_structure
    assert builder.vasp.code == mock_vasp
    assert builder.vasp.parameters is not None
    assert builder.relax_settings['algo']


@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_band_protocol(basic_env, mock_vasp, vasp_structure, potcar_family_name):
    """Test VASP band structure workchain protocol"""

    builder = VaspBandsWorkChain.get_builder_from_protocol(
        code=mock_vasp,
        structure=vasp_structure,
        overrides={
            'scf': {'potential_family': potcar_family_name, 'potential_mapping': {'In_d': 'In_d'}},
            'relax': {'vasp': {'potential_family': potcar_family_name, 'potential_mapping': {'In_d': 'In_d'}}},
        },
    )

    assert builder.structure == vasp_structure
    assert builder.scf.code == mock_vasp
    assert builder.scf.parameters is not None
    assert builder.band_settings is not None
    assert builder.relax.relax_settings is not None

    # No relax
    builder = VaspBandsWorkChain.get_builder_from_protocol(
        code=mock_vasp,
        structure=vasp_structure,
        run_relax=False,
        overrides={
            'scf': {'potential_family': potcar_family_name, 'potential_mapping': {'In_d': 'In_d'}},
            'relax': {'vasp': {'potential_family': potcar_family_name, 'potential_mapping': {'In_d': 'In_d'}}},
        },
    )

    assert builder.structure == vasp_structure
    assert builder.scf.code == mock_vasp
    assert builder.scf.parameters is not None
    assert builder.band_settings is not None
    assert not builder.relax.relax_settings

    builder = VaspHybridBandsWorkChain.get_builder_from_protocol(
        code=mock_vasp,
        structure=vasp_structure,
        overrides={
            'scf': {'potential_family': potcar_family_name, 'potential_mapping': {'In_d': 'In_d'}},
            'relax': {'vasp': {'potential_family': potcar_family_name, 'potential_mapping': {'In_d': 'In_d'}}},
        },
    )

    assert builder.structure == vasp_structure
    assert builder.scf.code == mock_vasp
    assert builder.scf.parameters is not None
    assert builder.band_settings is not None
    assert builder.relax.relax_settings is not None


@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_conv_protocol(basic_env, mock_vasp, vasp_structure, potcar_family_name):
    """Test VASP convergence test workchain protocol"""

    builder = VaspConvergenceWorkChain.get_builder_from_protocol(
        code=mock_vasp,
        structure=vasp_structure,
        overrides={'potential_family': potcar_family_name, 'potential_mapping': {'In_d': 'In_d'}},
    )

    assert builder.structure == vasp_structure
    assert builder.vasp.code == mock_vasp
    assert builder.vasp.parameters is not None
    assert builder.conv_settings is not None
