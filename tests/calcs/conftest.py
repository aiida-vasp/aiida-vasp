"""
Fixtures for the VASP calculations.

-----------------------------------
Here we set up different pytest fixtures that are used to represent various VASP
calculations on which one can for instance test parsing etc.
"""

import pathlib

import pytest
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.common.folders import SandboxFolder
from aiida.engine.utils import instantiate_process
from aiida.manage.manager import get_manager
from aiida.orm import Dict, List
from aiida.plugins import DataFactory

from aiida_vasp.calcs.base import VaspCalcBase
from aiida_vasp.calcs.neb import VaspNEBCalculation
from aiida_vasp.calcs.vasp import VaspCalculation
from aiida_vasp.calcs.vasp2w90 import Vasp2w90Calculation


@pytest.fixture()
def sandbox_folder():
    """Yield a `SandboxFolder` that can be used for tests where a Folder is needed."""

    with SandboxFolder() as folder:
        yield folder


@pytest.fixture()
def base_calc(aiida_profile, vasp_code):
    """An instance of a VaspCalcBase Process."""

    manager = get_manager()
    runner = manager.get_runner()
    inputs = AttributeDict()

    metadata = AttributeDict({'options': {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}})

    inputs.code = vasp_code
    inputs.metadata = metadata

    return instantiate_process(runner, VaspCalcBase, **inputs)


@pytest.fixture()
def vasp_calc(vasp_inputs):
    """An instance of a VaspCalculation Process."""

    def inner(inputs=None, settings=None):
        if inputs is None:
            inputs = vasp_inputs(settings)
        manager = get_manager()
        runner = manager.get_runner()

        return instantiate_process(runner, VaspCalculation, **inputs)

    return inner


@pytest.fixture()
def vasp_neb_inputs(aiida_profile, vasp_params, vasp_kpoints, vasp_structure, potentials, vasp_code):
    """Inputs dictionary for CalcJob Processes."""

    def inner(settings=None, parameters=None):
        inputs = AttributeDict()

        metadata = AttributeDict({'options': {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}})

        if settings is not None:
            inputs.settings = Dict(dict=settings)

        if isinstance(parameters, dict):
            parameters = orm.Dict(dict=parameters)

        if parameters is None:
            parameters = AttributeDict(vasp_params.get_dict())
            parameters['images'] = 3
            parameters = orm.Dict(dict=parameters)

        inputs.code = vasp_code
        inputs.metadata = metadata
        inputs.parameters = parameters
        inputs.kpoints, _ = vasp_kpoints

        inputs.initial_structure = vasp_structure
        inputs.final_structure = vasp_structure

        inputs.potential = potentials

        neb_images = {f'images_{idx:02d}': vasp_structure for idx in range(1, 4)}
        inputs.neb_images = neb_images

        return inputs

    return inner


@pytest.fixture()
def vasp_neb_calc(vasp_neb_inputs):
    """An instance of a VaspCalculation Process."""

    def inner(inputs=None, settings=None):
        if inputs is None:
            inputs = vasp_neb_inputs(settings)
        manager = get_manager()
        runner = manager.get_runner()

        return instantiate_process(runner, VaspNEBCalculation, **inputs)

    return inner


@pytest.fixture()
def vasp2w90_calc(vasp_inputs):
    """An instance of a VaspCalculation Process."""

    def inner(inputs=None, settings=None):
        if inputs is None:
            inputs = vasp_inputs(settings)

        manager = get_manager()
        runner = manager.get_runner()

        return instantiate_process(runner, Vasp2w90Calculation, **inputs)

    return inner


@pytest.fixture
def vasp_calc_and_ref(vasp_calc, vasp_kpoints, ref_incar):
    """Fixture for non varying setup of a vasp calculation."""
    calc = vasp_calc(settings={'parser_settings': {'add_bands': True, 'add_dos': True}})
    _, ref_kpoints = vasp_kpoints

    return calc, {'kpoints': ref_kpoints, 'incar': ref_incar}


@pytest.fixture
def vasp2w90_calc_and_ref(vasp2w90_calc, vasp_kpoints, vasp2w90_inputs, ref_incar_vasp2w90, ref_win):
    """Fixture for non varying setup of a vasp2w90 calculation."""

    inputs = vasp2w90_inputs(
        settings={
            'parser_settings': {
                'add_bands': True,
                'add_dos': True,
                'poscar_precision': 12,
            }
        }
    )

    calc = vasp2w90_calc(inputs=inputs)
    _, ref_kpoints = vasp_kpoints

    return calc, {'kpoints': ref_kpoints, 'incar': ref_incar_vasp2w90, 'win': ref_win}


@pytest.fixture()
def vasp_chgcar(aiida_profile, data_path):
    """CHGCAR node and reference fixture."""

    chgcar_path = data_path('chgcar', 'CHGCAR')
    chgcar = DataFactory('vasp.chargedensity')(file=chgcar_path)
    with open(chgcar_path, 'r', encoding='utf8') as ref_chgcar_fo:
        ref_chgcar = ref_chgcar_fo.read()
    return chgcar, ref_chgcar


@pytest.fixture()
def vasp_nscf_and_ref(vasp_calc_and_ref, vasp_chgcar, vasp_wavecar):
    """Fixture: vasp calc with chgcar and wavecar given."""
    calc, ref = vasp_calc_and_ref
    chgcar, ref_chgcar = vasp_chgcar
    wavecar, ref_wavecar = vasp_wavecar
    calc.use_charge_density(chgcar)
    calc.use_wavefunctions(wavecar)
    calc.inp.parameters.update_dict({'icharg': 11})
    ref['chgcar'] = ref_chgcar
    ref['wavecar'] = ref_wavecar
    return calc, ref


@pytest.fixture()
def vasp_inputs(aiida_profile, vasp_params, vasp_kpoints, vasp_structure, potentials, vasp_code):
    """Inputs dictionary for CalcJob Processes."""

    def inner(settings=None, parameters=None):
        inputs = AttributeDict()

        metadata = AttributeDict({'options': {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}})

        if settings is not None:
            inputs.settings = Dict(dict=settings)

        if isinstance(parameters, dict):
            parameters = orm.Dict(dict=parameters)

        if parameters is None:
            parameters = AttributeDict(vasp_params.get_dict())
            parameters = orm.Dict(dict=parameters)
        inputs.code = vasp_code
        inputs.metadata = metadata
        inputs.parameters = parameters
        inputs.kpoints, _ = vasp_kpoints
        inputs.structure = vasp_structure
        inputs.potential = potentials

        return inputs

    return inner


@pytest.fixture
def ref_incar(data_path):
    with open(data_path('incar', 'INCAR'), 'r', encoding='utf8') as reference_incar_fo:
        # yield reference_incar_fo.read().strip()
        yield reference_incar_fo.readlines()


@pytest.fixture()
def vasp_wavecar(aiida_profile, data_path):
    """WAVECAR node and reference fixture."""

    wavecar_path = data_path('wavecar', 'WAVECAR')
    wavecar = DataFactory('vasp.wavefun')(file=wavecar_path)
    with open(wavecar_path, 'r', encoding='utf8') as ref_wavecar_fo:
        ref_wavecar = ref_wavecar_fo.read()
    return wavecar, ref_wavecar


@pytest.fixture()
def vasp2w90_inputs(
    aiida_profile,
    vasp_params,
    vasp_kpoints,
    vasp_structure,
    potentials,
    vasp_code,
    wannier_projections,
    wannier_params,
):
    """Inputs dictionary for CalcJob Processes."""

    def inner(settings=None, parameters=None):
        inputs = AttributeDict()

        metadata = AttributeDict({'options': {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}})

        if settings is not None:
            inputs.settings = Dict(dict=settings)

        if isinstance(parameters, dict):
            parameters = orm.Dict(dict=parameters)

        if parameters is None:
            parameters = AttributeDict(vasp_params.get_dict())
            parameters = orm.Dict(dict=parameters)

        inputs.code = vasp_code
        inputs.metadata = metadata
        inputs.parameters = parameters
        inputs.kpoints, _ = vasp_kpoints
        inputs.structure = vasp_structure
        inputs.potential = potentials

        inputs.wannier_parameters = wannier_params
        inputs.wannier_projections = wannier_projections

        return inputs

    return inner


@pytest.fixture
def wannier_projections():
    wannier_projections = List()
    wannier_projections.extend(['Ga : s; px; py; pz', 'As : px; py; pz'])
    return wannier_projections


@pytest.fixture
def wannier_params():
    return Dict(
        dict=dict(  # pylint: disable=use-dict-literal
            dis_num_iter=1000,
            num_bands=24,
            num_iter=0,
            num_wann=14,
            spinors=True,
        )
    )


@pytest.fixture
def ref_incar_vasp2w90(data_path):
    with open(data_path('wannier', 'INCAR'), 'r', encoding='utf8') as reference_incar_wannier:
        yield reference_incar_wannier.readlines()


@pytest.fixture
def ref_win(data_path):
    data = pathlib.Path(data_path('wannier90.win'))
    yield data.read_text(encoding='utf8')
