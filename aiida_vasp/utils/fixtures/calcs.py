"""VaspCalculation fixtures"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import pytest

from aiida_vasp.utils.aiida_utils import get_data_node
from .data import vasp_code, vasp_params, potentials, vasp_kpoints, vasp_structure, ref_incar, vasp_chgcar, vasp_wavecar, wannier_params, \
    wannier_projections, ref_win


@pytest.fixture
def vasp_calc_and_ref(create_calc_and_ref, ref_incar):
    """Fixture for non varying setup of a vasp calculation"""
    from aiida_vasp.calcs.vasp import VaspCalculation
    calc, ref = create_calc_and_ref(VaspCalculation, ref_incar=ref_incar)
    calc.use_settings(get_data_node('parameter', dict={'parser_settings': {'add_bands': True, 'add_dos': True}}))
    return calc, ref


@pytest.fixture
def vasp2w90_calc_and_ref(create_calc_and_ref, ref_incar_vasp2w90, wannier_params, wannier_projections, ref_win):
    """Fixture for non varying setup of a vasp2w90 calculation"""
    from aiida_vasp.calcs.vasp2w90 import Vasp2w90Calculation
    calc, ref = create_calc_and_ref(Vasp2w90Calculation, ref_incar=ref_incar_vasp2w90)
    calc.use_wannier_parameters(wannier_params)
    calc.use_wannier_projections(wannier_projections)
    calc.use_settings(get_data_node('parameter', dict={'poscar_precision': 12}))
    ref['win'] = ref_win
    return calc, ref


@pytest.fixture
def create_calc_and_ref(vasp_code, vasp_params, potentials, vasp_kpoints, vasp_structure):
    """A function that creates a calculation of a given type."""

    def inner(calc_type, ref_incar):
        """Create a calculation of the given type."""
        calc = calc_type()
        calc.use_code(vasp_code)
        calc.set_computer(vasp_code.get_computer())
        calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        calc.use_parameters(vasp_params)
        calc.use_potential(potentials['In'], kind='In')
        calc.use_potential(potentials['As'], kind='As')
        calc.use_potential(potentials['In_d'], kind='In_d')
        calc.use_structure(vasp_structure)
        kpoints, ref_kpoints = vasp_kpoints
        calc.use_kpoints(kpoints)
        return calc, {'kpoints': ref_kpoints, 'incar': ref_incar}

    return inner


@pytest.fixture()
def vasp_nscf_and_ref(vasp_calc_and_ref, vasp_chgcar, vasp_wavecar):
    """Fixture: vasp calc with chgcar and wavecar given"""
    calc, ref = vasp_calc_and_ref
    chgcar, ref_chgcar = vasp_chgcar
    wavecar, ref_wavecar = vasp_wavecar
    calc.use_charge_density(chgcar)
    calc.use_wavefunctions(wavecar)
    calc.inp.parameters.update_dict({'icharg': 11})
    ref['chgcar'] = ref_chgcar
    ref['wavecar'] = ref_wavecar
    return calc, ref


ONLY_ONE_CALC = pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('cif', 'mesh')], indirect=True)

STRUCTURE_TYPES = pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('cif', 'mesh'), ('str', 'mesh')], indirect=True)
