# pylint: disable=attribute-defined-outside-init
"""
BandsWorkChain.

Intended to be used to extract the band structure using SeeKpath as a preprossesor
to extract the k-point path.
"""
import enum
from aiida.common.extendeddicts import AttributeDict
from aiida.engine.workchain import WorkChain, append_
from aiida.orm import WorkflowFactory
from aiida.engine.workfunctions import workfunction
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.utils.workchains import prepare_process_inputs


class BandsWorkChain(WorkChain):
    """Extract the band structure using k-point paths fetched from SeeKpath."""

    _verbose = False
    _next_workchain_string = 'vasp.vasp'
    _next_workchain = WorkflowFactory(_next_workchain_string)

    class ChargeEnum(enum.IntEnum):
        """Encode values for the initial charge density in enum."""
        WAVE = 0
        CHARGE = 1
        ATOMIC = 2
        CONSTANT_CHARGE = 11
        CONSTANT_ATOMIC = 12

    class OrbitEnum(enum.IntEnum):
        """Encode values for the projector information in enum."""
        ATOM = 0
        ATOM_LM = 1
        ATOM_LM_PHASE = 2
        NO_RWIGS_ATOM = 10
        NO_RWIGS_ATOM_LM = 11
        NO_RWIGS_ATOM_LM_PHASE = 12
        ATOM_LM_WAVE = 5

        @classmethod
        def get_from_combination(cls, **kwargs):
            """Get the correct mode of the projectors/decomposition."""
            combination = tuple(kwargs[i] for i in ['lm', 'phase', 'wigner_seitz_radius'])
            value_from_combinations = {
                (True, True, True): cls.ATOM_LM_PHASE,
                (False, False, False): cls.NO_RWIGS_ATOM,
                (False, True, True): cls.ATOM_LM_PHASE,
                (True, True, False): cls.NO_RWIGS_ATOM_LM_PHASE,
                (True, False, True): cls.ATOM_LM,
                (True, False, False): cls.NO_RWIGS_ATOM_LM,
                (False, False, True): cls.ATOM,
                (False, True, False): cls.NO_RWIGS_ATOM_LM_PHASE
            }
            return value_from_combinations[combination]

    @classmethod
    def define(cls, spec):
        super(BandsWorkChain, cls).define(spec)
        spec.expose_inputs(cls._next_workchain, exclude=('kpoints', 'parameters', 'settings'))
        spec.input('chgcar', valid_type=get_data_class('vasp.chargedensity'))
        spec.input('parameters', valid_type=get_data_class('dict'), required=False)
        spec.input('bands_parameters', valid_type=get_data_class('dict'), required=False)
        spec.input('settings', valid_type=get_data_class('dict'), required=False)
        spec.input(
            'kpoints_distance',
            valid_type=get_data_class('float'),
            required=False,
            default=get_data_node('float', 0.05),
            help="""
            The distance between each k-point along each high-symmetry line.
            """)
        spec.input(
            'decompose_bands',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
            Decompose the band structure on each atom.
            """)
        spec.input(
            'decompose_wave',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
            Decompose the wave function.
            """)
        spec.input(
            'lm',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
            Further decompose the decomposition into l- and m-states.
            """)
        spec.input(
            'phase',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
            Further decompose the l- and m-state decomposition into phases.
            """)
        spec.input(
            'wigner_seitz_radius',
            valid_type=get_data_class('list'),
            required=False,
            default=get_data_node('list', list=[False]),
            help="""
            The Wigner-Seitz radius for each atom type in AA as a list. If set, the internal projectors are not utilzed.
            """)
        spec.outline(
            cls.initialize,
            cls.get_kpoints_path,
            cls.init_next_workchain,
            cls.run_next_workchain,
            cls.verify_next_workchain,
            cls.results,
            cls.finalize
        )  # yapf: disable

        spec.expose_outputs(cls._next_workchain)
        spec.output('output_parameters_seekpath', valid_type=get_data_class('dict'))
        spec.output('output_bands', valid_type=get_data_class('array.bands'))
        spec.output('output_kpoints', valid_type=get_data_class('array.kpoints'))
        spec.output('output_structure_primitive', valid_type=get_data_class('structure'))

    def initialize(self):
        """Initialize."""
        self._init_context()
        self._init_inputs()
        self._init_settings()

        return

    def _init_context(self):
        """Initialize context variables."""
        self.ctx.inputs = AttributeDict()

        return

    def _init_settings(self):
        """Initialize the settings."""
        # Make sure we parse the bands
        if 'settings' in self.inputs:
            settings = AttributeDict(self.inputs.settings.get_dict())
        else:
            settings = AttributeDict({'parser_settings': {}})
        dict_entry = {'add_bands': True}
        try:
            settings.parser_settings.update(dict_entry)
        except AttributeError:
            settings.parser_settings = dict_entry
        self.ctx.inputs.settings = settings

        return

    def _init_inputs(self):
        """Initialize inputs."""
        self.ctx.inputs.parameters = self._init_parameters()

        self.ctx.inputs.seekpath_parameters = get_data_node('dict', dict={'reference_distance': self.inputs.kpoints_distance.value})

        try:
            self._verbose = self.inputs.verbose.value
        except AttributeError:
            pass

    def _add_overrides(self, parameters):
        """Add parameters tag overrides, except the ones controlled by other inputs (for provenance)."""
        overrides = AttributeDict({k.lower(): v for k, v in self.inputs.bands_parameters.get_dict().items()})
        parameters.update(overrides)

    def _init_parameters(self):
        """Set parameters based on either supplied or forced entries."""
        parameters = AttributeDict()
        try:
            input_parameters = self.inputs.parameters.get_dict()
            check_parameters_bands_entries(input_parameters)
            parameters.update(input_parameters)
        except AttributeError:
            pass
        if 'bands_parameters' in self.inputs:
            # Add override parameters (if user force what to use)
            try:
                self._add_overrides(parameters)
            except ValueError as err:
                self._fail_compat(exception=err)
        else:
            # Add plugin controlled flags
            self._set_icharg(parameters)
            self._set_lorbit(parameters)
            self._set_wigner_seitz_radius(parameters)
            self._set_ismear(parameters)

        return parameters

    def _set_ismear(self, parameters):
        """Make sure we do not supply invalid integration methods when running explicit k-point grids."""
        try:
            ismear = int(parameters.ismear)
            if ismear < -1:
                self.report('The requested integration method is incompatible with explicit k-point grids. '
                            'Setting the integration method to the default (removing tag).')
                del parameters.ismear
        except AttributeError:
            pass

    def _set_icharg(self, parameters):
        """Set the flag to start from input charge density and keep it constant."""
        parameters.icharg = self.ChargeEnum.CONSTANT_CHARGE

    def _set_wigner_seitz_radius(self, parameters):
        """Set the Wigner Seitz radius that is used to project/decompose."""
        wigner_seitz_radius = self.inputs.wigner_seitz_radius.get_list()
        if wigner_seitz_radius[0]:
            parameters.rwigs = wigner_seitz_radius

    def _set_lorbit(self, parameters):
        """Set the flag that controls the projectors/decomposition onto orbitals."""
        if self.inputs.decompose_bands:
            if self.inputs.decompose_wave:
                # Issue a warning that one can only use either or
                raise ValueError('Only projections/decompositions on the bands or the ' 'wave function are allowed.')
            wigner_seitz_radius = False
            if self.inputs.wigner_seitz_radius.value[0]:
                wigner_seitz_radius = True
            parameters.lorbit = self.OrbitEnum.get_from_combination(
                lm=self.inputs.lm.value, phase=self.inputs.phase.value, wigner_seitz_radius=wigner_seitz_radius)
        else:
            if self.inputs.decompose_wave:
                parameters.lorbit = self.OrbitEnum.ATOM_LM_WAVE

    def init_next_workchain(self):
        """Initialize the next workchain."""
        try:
            self.ctx.inputs
        except AttributeError:
            raise ValueError('No input dictionary was defined in self.ctx.inputs')

        # Add exposed inputs
        self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))

        # Make sure we do not have any floating dict (convert to Dict)
        self.ctx.inputs = prepare_process_inputs(self.ctx.inputs)

    def run_next_workchain(self):
        """Run the next workchain."""
        inputs = self.ctx.inputs
        running = self.submit(self._next_workchain, **inputs)

        if hasattr(running, 'pid'):
            self.report('launching {}<{}> '.format(self._next_workchain.__name__, running.pid))
        else:
            # Aiida < 1.0
            self.report('launching {}<{}> '.format(self._next_workchain.__name__, running.pk))

        return self.to_context(workchains=append_(running))

    def get_kpoints_path(self):
        """
        Fetch the k-point path.

        Run SeeKpath to get the high symmetry lines of the given structure. This
        routine returns a new (potentially different to the input structure) primitive
        structure. It also returns the k-point path for this structure.
        """
        result = seekpath_structure_analysis(self.inputs.structure, self.ctx.inputs.seekpath_parameters)

        self.ctx.inputs.structure = result['primitive_structure']
        self.ctx.inputs.kpoints = result['explicit_kpoints']

        # Set the output nodes for the primitive structure and the k-point path
        self.out('output_structure_primitive', result['primitive_structure'])
        self.out('output_kpoints', result['explicit_kpoints'])
        self.out('output_parameters_seekpath', result['parameters'])

        return

    def verify_next_workchain(self):
        """Verify and inherit exit status from child workchains."""

        workchain = self.ctx.workchains[-1]
        # Adopt exit status from last child workchain (supposed to be
        # successfull)
        next_workchain_exit_status = workchain.exit_status
        if not next_workchain_exit_status:
            self.exit_status = 0
        else:
            self.exit_status = next_workchain_exit_status
            self.report('The child {}<{}> returned a non-zero exit status, {}<{}> '
                        'inherits exit status {}'.format(workchain.__class__.__name__, workchain.pk, self.__class__.__name__, self.pid,
                                                         next_workchain_exit_status))
        return

    def results(self):
        """Attach the remaining output results."""

        if not self.exit_status:
            workchain = self.ctx.workchains[-1]
            self.out_many(self.exposed_outputs(workchain, self._next_workchain))

        return

    def finalize(self):
        """Finalize the workchain."""
        return self.exit_status

    def _fail_compat(self, *args, **kwargs):
        """Method to handle general failures."""
        if hasattr(self, 'fail'):
            self.fail(*args, **kwargs)  # pylint: disable=no-member
        else:
            msg = '{}'.format(kwargs['exception'])
            self.abort_nowait(msg)  # pylint: disable=no-member
        return


@workfunction
def seekpath_structure_analysis(structure, parameters):
    """
    Workfunction to extract k-points in the reciprocal cell.

    This workfunction will take a structure and pass it through SeeKpath to get the
    primitive cell and the path of high symmetry k-points through its Brillouin zone.
    Note that the returned primitive cell may differ from the original structure in
    which case the k-points are only congruent with the primitive cell.
    """
    from aiida.tools import get_explicit_kpoints_path
    return get_explicit_kpoints_path(structure, **parameters.get_dict())


def check_parameters_bands_entries(parameters):
    """Check that some flags are not present in the parameters (no override is allowed)."""

    overrides = AttributeDict({k.lower(): v for k, v in parameters.items()})
    if 'icharg' in overrides:
        raise ValueError('overriding ICHARG not allowed, use inputs to control')
    if 'lorbit' in overrides:
        raise ValueError('overriding ISIF not allowed, use inputs to control')
