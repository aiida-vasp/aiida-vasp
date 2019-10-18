""" # noqa: D205
Bands workchain
---------------
Intended to be used to extract the band structure using SeeKpath as a preprossesor
to extract the k-point path.
"""

# pylint: disable=attribute-defined-outside-init
import enum
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import WorkChain, append_, calcfunction
from aiida.plugins import WorkflowFactory
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.utils.workchains import prepare_process_inputs, compose_exit_code


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
        spec.expose_inputs(cls._next_workchain, exclude=('parameters', 'settings'))
        spec.input('parameters', valid_type=get_data_class('dict'), required=False)
        spec.input('bands.parameters', valid_type=get_data_class('dict'), required=False)
        spec.input('settings', valid_type=get_data_class('dict'), required=False)
        spec.input('bands.kpoints_distance',
                   valid_type=get_data_class('float'),
                   required=False,
                   default=get_data_node('float', 0.05),
                   help="""
            The distance between each k-point along each high-symmetry line.
            """)
        spec.input('bands.decompose_bands',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=get_data_node('bool', False),
                   help="""
            Decompose the band structure on each atom.
            """)
        spec.input('bands.decompose_wave',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=get_data_node('bool', False),
                   help="""
            Decompose the wave function.
            """)
        spec.input('bands.lm',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=get_data_node('bool', False),
                   help="""
            Further decompose the decomposition into l- and m-states.
            """)
        spec.input('bands.phase',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=get_data_node('bool', False),
                   help="""
            Further decompose the l- and m-state decomposition into phases.
            """)
        spec.input('bands.wigner_seitz_radius',
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
        spec.output('seekparam', valid_type=get_data_class('dict'))
        spec.output('bands', valid_type=get_data_class('array.bands'))
        spec.output('kpoints', valid_type=get_data_class('array.kpoints'))
        spec.output('structure_primitive', valid_type=get_data_class('structure'))
        spec.exit_code(0, 'NO_ERROR', message='the sun is shining')
        spec.exit_code(420, 'ERROR_NO_CALLED_WORKCHAIN', message='no called workchain detected')
        spec.exit_code(500, 'ERROR_UNKNOWN', message='unknown error detected in the bands workchain')

    def initialize(self):
        """Initialize."""
        self._init_context()
        self._init_inputs()
        self._init_settings()

    def _init_context(self):
        """Initialize context variables."""
        self.ctx.exit_code = self.exit_codes.ERROR_UNKNOWN  # pylint: disable=no-member
        self.ctx.inputs = AttributeDict()

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

    def _init_inputs(self):
        """Initialize inputs."""
        self.ctx.inputs.parameters = self._init_parameters()

        # Do not put the SeeKPath parameters in the inputs to avoid port checking
        # of the next workchain
        self.ctx.seekpath_parameters = get_data_node('dict', dict={'reference_distance': self.inputs.bands.kpoints_distance.value})

        try:
            self._verbose = self.inputs.verbose.value
        except AttributeError:
            pass

    def _add_overrides(self, parameters):
        """Add parameters tag overrides, except the ones controlled by other inputs (for provenance)."""
        overrides = AttributeDict({k.lower(): v for k, v in self.inputs.bands.parameters.get_dict().items()})
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
        wigner_seitz_radius = self.inputs.bands.wigner_seitz_radius.get_list()
        if wigner_seitz_radius[0]:
            parameters.rwigs = wigner_seitz_radius

    def _set_lorbit(self, parameters):
        """Set the flag that controls the projectors/decomposition onto orbitals."""
        if self.inputs.bands.decompose_bands:
            if self.inputs.bands.decompose_wave:
                # Issue a warning that one can only use either or
                raise ValueError('Only projections/decompositions on the bands or the ' 'wave function are allowed.')
            wigner_seitz_radius = False
            if self.inputs.bands.wigner_seitz_radius.value[0]:
                wigner_seitz_radius = True
            parameters.lorbit = self.OrbitEnum.get_from_combination(lm=self.inputs.bands.lm.value,
                                                                    phase=self.inputs.bands.phase.value,
                                                                    wigner_seitz_radius=wigner_seitz_radius)
        else:
            if self.inputs.bands.decompose_wave:
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
        result = seekpath_structure_analysis(self.inputs.structure, self.ctx.seekpath_parameters)

        self.ctx.inputs.structure = result['primitive_structure']
        self.report('explicit points:{}'.format(result['explicit_kpoints']))
        self.ctx.inputs.kpoints = result['explicit_kpoints']
        self.report('input kpoints:{}'.format(self.ctx.inputs.kpoints))

        # Set the output nodes for the primitive structure and the seekpath parameters
        #if self._verbose:
        #    self.report("attaching the node {}<{}> as '{}'".format(node.__class__.__name__, node.pk, name))  # pylint: disable=not-callable
        self.out('structure_primitive', result['primitive_structure'])
        self.out('seekparam', result['parameters'])

    def verify_next_workchain(self):
        """Verify and inherit exit status from child workchains."""

        try:
            workchain = self.ctx.workchains[-1]
        except IndexError:
            self.report('There is no {} in the called workchain list.'.format(self._next_workchain.__name__))
            return self.exit_codes.ERROR_NO_CALLED_WORKCHAIN  # pylint: disable=no-member

        # Inherit exit status from last workchain (supposed to be
        # successfull)
        next_workchain_exit_status = workchain.exit_status
        next_workchain_exit_message = workchain.exit_message
        if not next_workchain_exit_status:
            self.ctx.exit_code = self.exit_codes.NO_ERROR  # pylint: disable=no-member
        else:
            self.ctx.exit_code = compose_exit_code(next_workchain_exit_status, next_workchain_exit_message)
            self.report('The called {}<{}> returned a non-zero exit status. '
                        'The exit status {} is inherited'.format(workchain.__class__.__name__, workchain.pk, self.ctx.exit_code))

        return self.ctx.exit_code

    def results(self):
        """Attach the remaining output results."""

        workchain = self.ctx.workchains[-1]
        self.out_many(self.exposed_outputs(workchain, self._next_workchain))

    def finalize(self):
        """Finalize the workchain."""
        return self.ctx.exit_code


@calcfunction
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
