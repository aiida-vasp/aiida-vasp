"""VaspCalculation fixtures"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import pytest

from .data import vasp_code, vasp_params, potentials, vasp_kpoints, vasp_structure, ref_incar, vasp_chgcar, vasp_wavecar


@pytest.fixture()
def vasp_calc_and_ref(vasp_code, vasp_params, potentials, vasp_kpoints, vasp_structure, ref_incar):
    """Fixture for non varying setup of a vasp calculation"""
    from aiida_vasp.calcs.vasp import VaspCalculation
    calc = VaspCalculation()
    calc.use_code(vasp_code)
    calc.set_computer(vasp_code.get_computer())
    calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine': 1})
    calc.use_parameters(vasp_params)
    calc.use_potential(potentials['In'], kind='In')
    calc.use_potential(potentials['As'], kind='As')
    calc.use_structure(vasp_structure)
    kpoints, ref_kpoints = vasp_kpoints
    calc.use_kpoints(kpoints)
    return calc, {'kpoints': ref_kpoints, 'incar': ref_incar}


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
