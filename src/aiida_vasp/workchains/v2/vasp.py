"""
VASP workchain.

---------------
Contains the VaspWorkChain class definition which uses the BaseRestartWorkChain.
"""
import numpy as np

from aiida.common.exceptions import InputValidationError, NotExistent
# from aiida.engine.job_processes import override
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import if_, while_
from aiida.engine.processes.workchains.restart import BaseRestartWorkChain, process_handler
from aiida.orm import Code, Dict, KpointsData
from aiida.orm.nodes.data.base import to_aiida_type
from aiida.plugins import CalculationFactory

from aiida_vasp.assistant.parameters import ParametersMassage
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.utils.workchains import compose_exit_code, prepare_process_inputs
from aiida_vasp.workchains.vasp import VaspWorkChain as VanillaVaspWorkChain

from .common import parameters_validator, site_magnetization_to_magmom
from .inputset.vaspsets import get_ldau_keys

assert issubclass(
    VanillaVaspWorkChain, BaseRestartWorkChain
), 'vasp.vasp workchain is not a subclass of BaseRestartWorkChain from aiida-core'


class VaspWorkChain(VanillaVaspWorkChain):
    """
    The VASP workchain.

    -------------------
    Error handling enriched wrapper around VaspCalculation.

    Deliberately conserves most of the interface (required inputs) of the VaspCalculation class, but
    makes it possible for a user to interact with a workchain and not a calculation.

    This is intended to be used instead of directly submitting a VaspCalculation,
    so that future features like
    automatic restarting, error checking etc. can be propagated to higher level workchains
    automatically by implementing them here.

    Usage::

        from aiida.common.extendeddicts import AttributeDict
        from aiida.work import submit
        basevasp = WorkflowFactory('vasp.vasp')
        inputs = basevasp.get_builder()
        inputs = AttributeDict()
        ## ... set inputs
        submit(basevasp, **inputs)

    To see a working example, including generation of input nodes from scratch, please
    refer to ``examples/run_vasp_lean.py``.


    Additional functionalities:

    - Automatic setting LDA+U key using the ``ldau_mapping`` input port.

    - Set kpoints using spacing in A^-1 * 2pi with the ``kpoints_spacing`` input port.

    - Perform dryrun and set parameters such as KPAR and NCORE automatically if ``auto_parallel`` input port exists.
      this will give rise to an additional output node ``parallel_settings`` containing the strategy obtained.

    """

    _verbose = False
    _calculation = CalculationFactory('vasp.vasp')

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.input('code', valid_type=Code)
        spec.input(
            'structure',
            valid_type=(get_data_class('structure'), get_data_class('cif')),
            required=True,
        )
        spec.input('kpoints', valid_type=get_data_class('array.kpoints'), required=False)
        spec.input(
            'potential_family',
            valid_type=get_data_class('str'),
            required=True,
            serializer=to_aiida_type,
        )
        spec.input(
            'potential_mapping',
            valid_type=get_data_class('dict'),
            required=True,
            serializer=to_aiida_type,
        )
        spec.input(
            'parameters',
            valid_type=get_data_class('dict'),
            required=True,
            validator=parameters_validator,
        )
        spec.input(
            'options',
            valid_type=get_data_class('dict'),
            required=True,
            serializer=to_aiida_type,
        )
        spec.input(
            'settings',
            valid_type=get_data_class('dict'),
            required=False,
            serializer=to_aiida_type,
        )
        spec.input('wavecar', valid_type=get_data_class('vasp.wavefun'), required=False)
        spec.input('chgcar', valid_type=get_data_class('vasp.chargedensity'), required=False)
        spec.input(
            'site_magnetization',
            valid_type=get_data_class('dict'),
            required=False,
            help='Site magnetization to be used as MAGMOM',
        )
        spec.input(
            'restart_folder',
            valid_type=get_data_class('remote'),
            required=False,
            help="""
            The restart folder from a previous workchain run that is going to be used.
            """,
        )
        spec.input(
            'max_iterations',
            valid_type=get_data_class('int'),
            required=False,
            default=lambda: get_data_node('int', 5),
            serializer=to_aiida_type,
            help="""
            The maximum number of iterations to perform.
            """,
        )
        spec.input(
            'clean_workdir',
            valid_type=get_data_class('bool'),
            required=False,
            serializer=to_aiida_type,
            default=lambda: get_data_node('bool', True),
            help="""
            If True, clean the work dir upon the completion of a successfull calculation.
            """,
        )
        spec.input(
            'verbose',
            valid_type=get_data_class('bool'),
            required=False,
            serializer=to_aiida_type,
            default=lambda: get_data_node('bool', False),
            help="""
            If True, enable more detailed output during workchain execution.
            """,
        )
        spec.input(
            'ldau_mapping',
            valid_type=get_data_class('dict'),
            required=False,
            serializer=to_aiida_type,
            help="""Settings for assign LDA+U related settings according to the input structure.

    mapping: a dictionary in the format of  {"Mn": [d, 4]...} for U
    utype: the type of LDA+U, default to 2, which is the one with only one parameter
    jmapping: a dictionary in the format of  {"Mn": [d, 4]...} but for J
    felec: Wether we are dealing with f electrons, will increase lmaxmix if we are.""",
        )
        spec.input(
            'kpoints_spacing',
            valid_type=get_data_class('float'),
            required=False,
            serializer=to_aiida_type,
            help='Spacing for the kpoints in units A^-1 * 2pi',
        )
        spec.input(
            'auto_parallel',
            valid_type=get_data_class('dict'),
            serializer=to_aiida_type,
            required=False,
            help='Automatic parallelisation settings, keywords passed to `get_jobscheme` function.',
        )
        spec.input(
            'dynamics.positions_dof',
            valid_type=get_data_class('list'),
            serializer=to_aiida_type,
            required=False,
            help="""
            Site dependent flag for selective dynamics when performing relaxation
            """,
        )
        spec.outline(
            cls.setup,
            cls.init_inputs,
            if_(cls.run_auto_parallel)(cls.prepare_inputs, cls.perform_autoparallel),
            while_(cls.should_run_process)(
                cls.prepare_inputs,
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )  # yapf: disable
        spec.output('parallel_settings', valid_type=get_data_class('dict'), required=False)

    def init_inputs(self):
        """Make sure all the required inputs are there and valid, create input dictionary for calculation."""

        #### START OF THE COPY FROM VASPWorkChain ####
        #  - the only change is that the section about kpoints is deleted
        self.ctx.inputs = AttributeDict()
        self.ctx.inputs.parameters = self._init_parameters()
        # Set the code
        self.ctx.inputs.code = self.inputs.code

        # Set the structure (poscar)
        self.ctx.inputs.structure = self.inputs.structure

        # Set the kpoints (kpoints) - No longer needed, using kpoints/kpoints spacing as from below
        # self.ctx.inputs.kpoints = self.inputs.kpoints

        # Set settings
        unsupported_parameters = None
        skip_parameters_validation = False
        if self.inputs.get('settings'):
            self.ctx.inputs.settings = self.inputs.settings
            # Also check if the user supplied additional tags that is not in the supported file.
            settings_dict = self.ctx.inputs.settings.get_dict()
            unsupported_parameters = settings_dict.get('unsupported_parameters', unsupported_parameters)
            skip_parameters_validation = settings_dict.get('skip_parameters_validation', skip_parameters_validation)

        # Perform inputs massage to accommodate generalization in higher lying workchains
        # and set parameters.
        try:
            parameters_massager = ParametersMassage(
                self.ctx.inputs.parameters,
                unsupported_parameters,
                skip_parameters_validation=skip_parameters_validation,
            )
        except Exception as exception:  # pylint: disable=broad-except
            return self.exit_codes.ERROR_IN_PARAMETER_MASSAGER.format(exception=exception)  # pylint: disable=no-member
        try:
            # Only set if they exists
            # Set any INCAR tags
            self.ctx.inputs.parameters = parameters_massager.parameters.incar
            # Set any dynamics input (currently only for selective dynamics, e.g. custom write to POSCAR)
            self.ctx.inputs.dynamics = parameters_massager.parameters.dynamics
            # Here we could set additional override flags, but those are not relevant for this VASP plugin
        except AttributeError:
            pass

        # Set options
        # Options is very special, not storable and should be
        # wrapped in the metadata dictionary, which is also not storable
        # and should contain an entry for options
        if 'options' in self.inputs:
            options = {}
            options.update(self.inputs.options)
            self.ctx.inputs.metadata = {'options': options}
            # Override the parser name if it is supplied by the user.
            parser_name = self.ctx.inputs.metadata['options'].get('parser_name')
            if parser_name:
                self.ctx.inputs.metadata['options']['parser_name'] = parser_name
            # Set MPI to True, unless the user specifies otherwise
            withmpi = self.ctx.inputs.metadata['options'].get('withmpi', True)
            self.ctx.inputs.metadata['options']['withmpi'] = withmpi

        # Utilise default input/output selections
        self.ctx.inputs.metadata['options']['input_filename'] = 'INCAR'
        self.ctx.inputs.metadata['options']['output_filename'] = 'OUTCAR'

        # Make sure we also bring along any label and description set on the WorkChain to the CalcJob, it if does
        # not exists, set to empty string.
        if 'metadata' in self.inputs:
            label = self.inputs.metadata.get('label', '')
            description = self.inputs.metadata.get('description', '')
            if 'metadata' not in self.ctx.inputs:
                self.ctx.inputs.metadata = {}
            self.ctx.inputs.metadata['label'] = label
            self.ctx.inputs.metadata['description'] = description

        # Carry on site magnetization for initialization
        if 'site_magnetization' in self.inputs:
            magmom = site_magnetization_to_magmom(self.inputs.site_magnetization.get_dict())
            assert len(magmom) == len(self.inputs.structure.sites)
            self.ctx.inputs.parameters['magmom'] = magmom

        # Verify and set potentials (potcar)
        if not self.inputs.potential_family.value:
            self.report('An empty string for the potential family name was detected.')  # pylint: disable=not-callable
            return self.exit_codes.ERROR_NO_POTENTIAL_FAMILY_NAME  # pylint: disable=no-member
        try:
            self.ctx.inputs.potential = get_data_class('vasp.potcar').get_potcars_from_structure(
                structure=self.inputs.structure,
                family_name=self.inputs.potential_family.value,
                mapping=self.inputs.potential_mapping.get_dict(),
            )
        except ValueError as err:
            return compose_exit_code(self.exit_codes.ERROR_POTENTIAL_VALUE_ERROR.status, str(err))  # pylint: disable=no-member
        except NotExistent as err:
            return compose_exit_code(self.exit_codes.ERROR_POTENTIAL_DO_NOT_EXIST.status, str(err))  # pylint: disable=no-member

        # Store verbose parameter in ctx - otherwise it will not work after deserialization
        try:
            self.ctx.verbose = self.inputs.verbose.value
        except AttributeError:
            self.ctx.verbose = self._verbose

        # Set the charge density (chgcar)
        if 'chgcar' in self.inputs:
            self.ctx.inputs.charge_density = self.inputs.chgcar

        # Set the wave functions (wavecar)
        if 'wavecar' in self.inputs:
            self.ctx.inputs.wavefunctions = self.inputs.wavecar

        ##### END OF THE COPY from VaspWorkChain   #####

        # Set the kpoints (kpoints)
        if 'kpoints' in self.inputs:
            self.ctx.inputs.kpoints = self.inputs.kpoints
        elif 'kpoints_spacing' in self.inputs:
            kpoints = KpointsData()
            kpoints.set_cell_from_structure(self.ctx.inputs.structure)
            kpoints.set_kpoints_mesh_from_density(self.inputs.kpoints_spacing.value * np.pi * 2)
            self.ctx.inputs.kpoints = kpoints
        else:
            raise InputValidationError("Must supply either 'kpoints' or 'kpoints_spacing'")

        # Setup LDAU keys
        if 'ldau_mapping' in self.inputs:
            ldau_settings = self.inputs.ldau_mapping.get_dict()
            ldau_keys = get_ldau_keys(self.ctx.inputs.structure, **ldau_settings)
            # Directly update the raw inputs passed to VaspCalculation
            self.ctx.inputs.parameters.update(ldau_keys)
        return None

    def run_auto_parallel(self):
        """Wether we should run auto-parallelisation test"""
        return ('auto_parallel' in self.inputs and self.inputs.auto_parallel.value is True)

    def perform_autoparallel(self):
        """Dry run and obtain the best parallelisation settings"""
        from .common.dryrun import get_jobscheme

        self.report('Performing local dryrun for auto-parallelisation')  # pylint: disable=not-callable

        ind = prepare_process_inputs(self.ctx.inputs)

        nprocs = self.ctx.inputs.metadata['options']['resources']['tot_num_mpiprocs']

        # Take the settings pass it to the function
        kwargs = self.inputs.auto_parallel.get_dict()
        if 'cpus_per_node' not in kwargs:
            kwargs['cpus_per_node'] = self.inputs.code.computer.get_default_mpiprocs_per_machine()

        # If the dryrun errored, proceed the workchain
        try:
            scheme = get_jobscheme(ind, nprocs, **kwargs)
        except Exception as error:
            self.report(f"Dry-run errorred, process with cautions, message: {error.args}")  # pylint: disable=not-callable
            return

        if (scheme.ncore is None) or (scheme.kpar is None):
            self.report(f"Error NCORE: {scheme.ncore}, KPAR: {scheme.kpar}")  # pylint: disable=not-callable
            return

        parallel_opts = {'ncore': scheme.ncore, 'kpar': scheme.kpar}
        self.report(f"Found optimum KPAR={scheme.kpar}, NCORE={scheme.ncore}")  # pylint: disable=not-callable
        self.ctx.inputs.parameters.update(parallel_opts)
        self.out(
            'parallel_settings',
            Dict(dict={
                'ncore': scheme.ncore,
                'kpar': scheme.kpar
            }).store(),
        )

    # In this workchain variant we default to ignore the NELM breaches in the middle of the calculation
    @process_handler(priority=850, enabled=True)
    def ignore_nelm_breach_relax(self, node):
        """
        Not a actual handler but works as a switch to bypass checks for NELM breaches in the middle of an ionic relaxation.
        """
        _ = node
        self.ctx.ignore_transient_nelm_breach = True
