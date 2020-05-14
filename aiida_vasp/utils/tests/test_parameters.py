"""Test aiida_parameters."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import,no-member

import pytest

from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.parameters import ParametersMassage


@pytest.fixture
def init_relax_parameters():
    """Fixture for a set of general input parameters for relaxation."""
    general_parameters = AttributeDict()
    general_parameters.relax = AttributeDict()
    general_parameters.relax.algo = 'cg'
    general_parameters.relax.energy_cutoff = 0.01
    general_parameters.relax.force_cutoff = 0.01
    general_parameters.relax.steps = 60
    general_parameters.relax.positions = True
    general_parameters.relax.shape = True
    general_parameters.relax.volume = True

    return general_parameters


@pytest.fixture
def init_smearing_parameters():
    """Fixture for a set of general input parameters regarding Fermi integral smearing."""
    general_parameters = AttributeDict()
    general_parameters.smearing = AttributeDict()
    general_parameters.smearing.gaussian = True

    return general_parameters


@pytest.fixture
def init_charge_parameters():
    """Fixture for a set of general input parameters regarding how to construct initial charge density."""
    general_parameters = AttributeDict()
    general_parameters.charge = AttributeDict()
    general_parameters.charge.constant_charge = True

    return general_parameters


@pytest.fixture
def init_bands_parameters():
    """Fixture for a set of general input parameters for band structure calculations."""
    general_parameters = AttributeDict()
    general_parameters.bands = AttributeDict()
    general_parameters.bands.decompose_bands = False
    general_parameters.bands.decompose_wave = False

    return general_parameters


@pytest.fixture
def init_simple_workchain():
    """Fixture to simulate a fake workchain to store the exit codes and a dummy report function."""
    workchain = AttributeDict()
    workchain.exit_codes = AttributeDict()
    workchain.exit_codes.ERROR_INVALID_PARAMETER_DETECTED = 1
    workchain.exit_codes.ERROR_MISSING_PARAMETER_DETECTED = 1
    workchain.report = print

    return workchain


def test_relax_parameters_all_set(init_relax_parameters):
    """Test all standard relaxation parameters are set."""
    massager = ParametersMassage(None, init_relax_parameters)
    assert massager.exit_code is None
    parameters = massager.parameters
    assert parameters.ediffg == -0.01
    assert parameters.ibrion == 2
    assert parameters.nsw == 60
    assert parameters.isif == 3


def test_relax_parameters_energy(init_relax_parameters):
    """Test no provided force cutoff. It should be set to a default value."""
    del init_relax_parameters.relax.force_cutoff
    massager = ParametersMassage(None, init_relax_parameters)
    assert massager.exit_code is None
    parameters = massager.parameters
    assert parameters.ediffg == 0.01


def test_relax_parameters_no_algo(init_relax_parameters, init_simple_workchain):
    """Test no provided algo tag."""
    mock_workchain = init_simple_workchain
    del init_relax_parameters.relax.algo
    massager = ParametersMassage(mock_workchain, init_relax_parameters)
    assert massager.exit_code is not None


def test_relax_parameters_vol_shape(init_relax_parameters):
    """Test volume and shape relaxation combinations."""
    del init_relax_parameters.relax.positions
    massager = ParametersMassage(None, init_relax_parameters)
    assert massager.exit_code is None
    parameters = massager.parameters
    assert parameters.isif == 6


def test_relax_parameters_pos_shape(init_relax_parameters):
    """Test position and shape relxation combinations."""
    del init_relax_parameters.relax.volume
    massager = ParametersMassage(None, init_relax_parameters)
    assert massager.exit_code is None
    parameters = massager.parameters
    assert parameters.isif == 4


def test_relax_parameters_vol(init_relax_parameters):
    """Test only volume relaxation."""
    del init_relax_parameters.relax.positions
    del init_relax_parameters.relax.shape
    massager = ParametersMassage(None, init_relax_parameters)
    assert massager.exit_code is None
    parameters = massager.parameters
    assert parameters.isif == 7


def test_relax_parameters_pos(init_relax_parameters):
    """Test only position relaxation."""
    del init_relax_parameters.relax.volume
    del init_relax_parameters.relax.shape
    massager = ParametersMassage(None, init_relax_parameters)
    assert massager.exit_code is None
    parameters = massager.parameters
    assert parameters.isif == 2


def test_relax_parameters_shape(init_relax_parameters):
    """Test only shape relaxation."""
    del init_relax_parameters.relax.volume
    del init_relax_parameters.relax.positions
    massager = ParametersMassage(None, init_relax_parameters)
    assert massager.exit_code is None
    parameters = massager.parameters
    assert parameters.isif == 5


def test_relax_parameters_nothing(init_relax_parameters):
    """Test if no relaxation parameters for volume, positions and shape are given."""
    del init_relax_parameters.relax.volume
    del init_relax_parameters.relax.positions
    del init_relax_parameters.relax.shape
    massager = ParametersMassage(None, init_relax_parameters)
    assert massager.exit_code is None
    parameters = massager.parameters
    assert parameters == AttributeDict()


def test_relax_parameters_override(init_relax_parameters):
    """Test what happens if we override a parameters."""
    value = 1
    init_relax_parameters.isif = value
    massager = ParametersMassage(None, init_relax_parameters)
    assert massager.exit_code is None
    parameters = massager.parameters
    assert parameters.isif == value


def test_smearing_parameters():
    """Test smearing parameters."""
    parameters = AttributeDict()
    parameters.smearing = AttributeDict()
    parameters.smearing.gaussian = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.ismear == 0
    parameters.smearing.gaussian = False
    parameters.smearing.fermi = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.ismear == -1
    parameters.smearing.fermi = False
    parameters.smearing.tetra = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.ismear == -5
    parameters.smearing.tetra = False
    parameters.smearing.mp = 4
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.ismear == 4


def test_charge_parameters():
    """Test charge parameters."""
    parameters = AttributeDict()
    parameters.charge = AttributeDict()
    parameters.charge.from_wave = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    parameters.charge.from_wave = False
    parameters.charge.from_charge = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.icharg == 1
    parameters.charge.from_charge = False
    parameters.charge.from_atomic = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.icharg == 2
    parameters.charge.from_atomic = False
    parameters.charge.from_potential = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.icharg == 4
    parameters.charge.from_potential = False
    parameters.charge.constant_charge = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.icharg == 11
    parameters.charge.constant_charge = False
    parameters.charge.constant_atomic = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.icharg == 12


def test_orbital_projections():  # pylint: disable=too-many-statements
    """Test the parameters associated with orbital projections."""
    parameters = AttributeDict()
    parameters.bands = AttributeDict()
    parameters.bands.decompose_wave = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.lorbit == 5
    parameters.bands.decompose_wave = False
    parameters.bands.decompose_bands = True
    parameters.bands.decompose_auto = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.lorbit == 14
    parameters.bands.decompose_auto = False
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.lorbit == 10
    parameters.bands.lm = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.lorbit == 11
    parameters.bands.phase = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.lorbit == 12
    parameters.bands.lm = False
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.lorbit == 12

    # Now do the once with a Wigner-Seitz radius supplied
    parameters.bands.wigner_seitz_radius = [2.0]
    parameters.bands.lm = False
    parameters.bands.phase = False
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.lorbit == 0
    print(massager.parameters.rwigs)
    assert int(massager.parameters.rwigs[0]) == 2
    parameters.bands.lm = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.lorbit == 1
    parameters.bands.phase = True
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.lorbit == 2
    parameters.bands.lm = False
    massager = ParametersMassage(None, parameters)
    assert massager.exit_code is None
    assert massager.parameters.lorbit == 2

    # Should raise ValueError if Wigner-Seitz radius is not defined as a list.
    parameters.bands.wigner_seitz_radius = 2.0
    with pytest.raises(ValueError):
        massager = ParametersMassage(None, parameters)
