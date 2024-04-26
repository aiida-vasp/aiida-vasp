"""
Fixtures for the VASP calculations.

-----------------------------------
Here we set up different pytest fixtures that are used to represent various VASP
calculations on which one can for instance test parsing etc.
"""

# pylint: disable=unused-import,unused-argument,redefined-outer-name, import-outside-toplevel
import pytest
from aiida.common.extendeddicts import AttributeDict
from aiida.engine.utils import instantiate_process
from aiida.manage.manager import get_manager
from aiida.orm import Code

from aiida_vasp.utils.aiida_utils import create_authinfo, get_data_class, get_data_node

from .data import (
    POTCAR_FAMILY_NAME,
    POTCAR_MAP,
)


@pytest.fixture()
def sandbox_folder():
    """Yield a `SandboxFolder` that can be used for tests where a Folder is needed."""
    from aiida.common.folders import SandboxFolder

    with SandboxFolder() as folder:
        yield folder


@pytest.fixture()
def calc_with_retrieved(localhost):
    """A rigged CalcJobNode for testing the parser and that the calculation retrieve what is expected."""
    from aiida.common.links import LinkType
    from aiida.orm import CalcJobNode  # pylint: disable=no-name-in-module
    from aiida.plugins import DataFactory

    def _inner(file_path, input_settings=None):
        # Create a test computer
        computer = localhost

        process_type = 'aiida.calculations:vasp.vasp'

        node = CalcJobNode(computer=computer, process_type=process_type)
        node.base.attributes.set('input_filename', 'INCAR')
        node.base.attributes.set('output_filename', 'OUTCAR')
        # node.set_attribute('error_filename', 'aiida.err')
        node.base.attributes.set('scheduler_stderr', '_scheduler-stderr.txt')
        node.base.attributes.set('scheduler_stdout', '_scheduler-stdout.txt')
        node.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        node.set_option('max_wallclock_seconds', 1800)

        if input_settings is None:
            input_settings = {}

        settings = DataFactory('core.dict')(dict=input_settings)
        node.base.links.add_incoming(settings, link_type=LinkType.INPUT_CALC, link_label='settings')
        settings.store()
        node.store()

        # Create a `FolderData` that will represent the `retrieved` folder. Store the test
        # output fixture in there and link it.
        retrieved = DataFactory('core.folder')()
        retrieved.base.repository.put_object_from_tree(file_path)
        retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='retrieved')
        retrieved.store()

        return node

    return _inner


@pytest.fixture()
def neb_calc_with_retrieved(localhost):
    """A rigged CalcJobNode for testing the parser and that the calculation retrieve what is expected."""
    from aiida.common.links import LinkType
    from aiida.orm import CalcJobNode  # pylint: disable=no-name-in-module
    from aiida.plugins import DataFactory

    def _inner(file_path, input_settings=None, nimgs=3):
        # Create a test computer
        computer = localhost

        process_type = 'aiida.calculations:vasp.vasp'

        node = CalcJobNode(computer=computer, process_type=process_type)
        node.base.attributes.set('input_filename', 'INCAR')
        node.base.attributes.set('output_filename', 'OUTCAR')
        # node.set_attribute('error_filename', 'aiida.err')
        node.base.attributes.set('scheduler_stderr', '_scheduler-stderr.txt')
        node.base.attributes.set('scheduler_stdout', '_scheduler-stdout.txt')
        node.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        node.set_option('max_wallclock_seconds', 1800)

        if input_settings is None:
            input_settings = {}

        settings = DataFactory('core.dict')(dict=input_settings)
        node.base.links.add_incoming(settings, link_type=LinkType.INPUT_CALC, link_label='settings')
        settings.store()

        # Add inputs with the number of images
        param = get_data_class('core.dict')(dict={'images': nimgs})
        node.base.links.add_incoming(param, link_type=LinkType.INPUT_CALC, link_label='parameters')
        param.store()

        node.store()

        # Create a `FolderData` that will represent the `retrieved` folder. Store the test
        # output fixture in there and link it.
        retrieved = DataFactory('core.folder')()
        retrieved.base.repository.put_object_from_tree(file_path)
        retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='retrieved')
        retrieved.store()

        return node

    return _inner


@pytest.fixture()
def base_calc(fresh_aiida_env, vasp_code):
    """An instance of a VaspCalcBase Process."""
    from aiida_vasp.calcs.base import VaspCalcBase

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
    from aiida_vasp.calcs.vasp import VaspCalculation

    def inner(inputs=None, settings=None):
        if inputs is None:
            inputs = vasp_inputs(settings)
        manager = get_manager()
        runner = manager.get_runner()

        return instantiate_process(runner, VaspCalculation, **inputs)

    return inner


@pytest.fixture()
def vasp_neb_calc(vasp_neb_inputs):
    """An instance of a VaspCalculation Process."""
    from aiida_vasp.calcs.neb import VaspNEBCalculation

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
    from aiida_vasp.calcs.vasp2w90 import Vasp2w90Calculation

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
def run_vasp_process(fresh_aiida_env, vasp_params, potentials, vasp_kpoints, vasp_structure, mock_vasp):
    """Setup a standard VaspCalculation or VaspWorkChain with the mock executable that accepts input overrides."""

    def inner(inputs=None, settings=None, test_case=None, process_type='calcjob'):
        """
        Run a VaspCalculation or VaspWorkChain with specified input and settings overrides.

        Specific outputs can be selected using the test_case parameter.

        The type of process is set with the process_type parameter.
        """
        from aiida.engine import run

        inpts = AttributeDict()
        inpts.structure = vasp_structure
        parameters = vasp_params.get_dict()
        options = {
            'withmpi': False,
            'queue_name': 'None',
            'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1},
            'max_wallclock_seconds': 3600,
        }
        if test_case is not None:
            # Allow to fetch special tests cases using the mock-vasp executable
            parameters['system'] = f'test-case:{test_case}'
        if process_type == 'calcjob':
            from aiida.plugins import CalculationFactory

            process = CalculationFactory('vasp.vasp')
            inpts.potential = get_data_class('vasp.potcar').get_potcars_from_structure(
                structure=inpts.structure,
                family_name=POTCAR_FAMILY_NAME,
                mapping=POTCAR_MAP,
            )
            inpts.parameters = get_data_class('core.dict')(dict=parameters)
            inpts.metadata = {}
            inpts.metadata['options'] = options
        elif process_type == 'workchain':
            from aiida.plugins import WorkflowFactory

            process = WorkflowFactory('vasp.vasp')
            inpts.potential_family = get_data_node('core.str', POTCAR_FAMILY_NAME)
            inpts.potential_mapping = get_data_node('core.dict', dict=POTCAR_MAP)
            inpts.parameters = get_data_node('core.dict', dict={'incar': parameters})
            inpts.options = get_data_node('core.dict', dict=options)
            inpts.max_iterations = get_data_node('core.int', 1)
            inpts.clean_workdir = get_data_node('core.bool', False)
            inpts.verbose = get_data_node('core.bool', True)
        else:
            raise ValueError(
                f"The supplied process_type: {process_type} is not supported. Use either 'calcjob' or 'workchain.'"
            )

        mock_vasp.store()
        create_authinfo(computer=mock_vasp.computer, store=True)
        inpts.code = Code.get_from_string('mock-vasp@localhost')
        kpoints, _ = vasp_kpoints
        inpts.kpoints = kpoints
        if inputs is not None:
            # Allow overrides of the input
            inpts.update(inputs)
        if settings is not None and isinstance(settings, dict):
            inpts.settings = get_data_node('core.dict', dict=settings)
        results_and_node = run.get_node(process, **inpts)
        return results_and_node

    return inner


ONLY_ONE_CALC = pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('cif', 'mesh')], indirect=True)

STRUCTURE_TYPES = pytest.mark.parametrize(
    ['vasp_structure', 'vasp_kpoints'],
    [('cif', 'mesh'), ('str', 'mesh')],
    indirect=True,
)
