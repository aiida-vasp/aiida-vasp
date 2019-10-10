""" # noqa: D205
Relax workchain
---------------
Structure relaxation workchain.
"""
# pylint: disable=attribute-defined-outside-init
import enum
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import WorkChain, append_, while_, if_
from aiida.plugins import WorkflowFactory

from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.utils.workchains import compare_structures, prepare_process_inputs, compose_exit_code


class RelaxWorkChain(WorkChain):
    """Structure relaxation workchain."""
    _verbose = False
    _next_workchain_string = 'vasp.verify'
    _next_workchain = WorkflowFactory(_next_workchain_string)

    class AlgoEnum(enum.IntEnum):
        """Encode values for algorithm descriptively in enum."""
        NO_UPDATE = -1
        IONIC_RELAXATION_RMM_DIIS = 1
        IONIC_RELAXATION_CG = 2

    class ModeEnum(enum.IntEnum):
        """Encode values for mode of relaxation descriptively in enum."""
        POS_ONLY = 1
        POS_SHAPE_VOL = 3
        POS_SHAPE = 4
        SHAPE_ONLY = 5

        @classmethod
        def get_from_dof(cls, **kwargs):
            """Get the correct mode of relaxation for the given degrees of freedom."""
            dof = tuple(kwargs[i] for i in ['positions', 'shape', 'volume'])
            value_from_dof = {
                (True, False, False): cls.POS_ONLY,
                (True, True, True): cls.POS_SHAPE_VOL,
                (True, True, False): cls.POS_SHAPE,
                (False, True, False): cls.SHAPE_ONLY
            }
            return value_from_dof[dof]

    @classmethod
    def define(cls, spec):
        super(RelaxWorkChain, cls).define(spec)
        spec.expose_inputs(cls._next_workchain, exclude=('parameters', 'structure', 'settings'))
        spec.input('structure', valid_type=(get_data_class('structure'), get_data_class('cif')))
        spec.input('parameters', valid_type=get_data_class('dict'), required=False)
        spec.input('relax.parameters', valid_type=get_data_class('dict'), required=False)
        spec.input('settings', valid_type=get_data_class('dict'), required=False)
        spec.input(
            'relax.perform',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
            If True, perform relaxation.
            """)
        spec.input(
            'relax.energy_cutoff',
            valid_type=get_data_class('float'),
            required=False,
            help="""
            The cutoff that determines when the relaxation procedure is stopped. In this
            case it stops when the total energy between two ionic steps is less than the
            supplied value.
            """)
        spec.input(
            'relax.force_cutoff',
            valid_type=get_data_class('float'),
            required=False,
            help="""
            The cutoff that determines when the relaxation procedure is stopped. In this
            case it stops when all forces are smaller than than the
            supplied value.
            """)
        spec.input(
            'relax.steps',
            valid_type=get_data_class('int'),
            required=False,
            default=get_data_node('int', 60),
            help="""
                   The number of relaxation steps to perform (updates to the atomic positions,
            unit cell size or shape).
                   """)
        spec.input(
            'relax.positions',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', True),
            help="""
                   If True, perform relaxation of the atomic positions.
                   """)
        spec.input(
            'relax.shape',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
                   If True, perform relaxation of the unit cell shape.
                   """)
        spec.input(
            'relax.volume',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
                   If True, perform relaxation of the unit cell volume..
                   """)
        spec.input(
            'relax.convergence_on',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
                   If True, test convergence based on selected criterias set.
                   """)
        spec.input(
            'relax.convergence_absolute',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
                   If True, test convergence based on absolute differences.
                   """)
        spec.input(
            'relax.convergence_max_iterations',
            valid_type=get_data_class('int'),
            required=False,
            default=get_data_node('int', 5),
            help="""
                   The number of iterations to perform if the convergence criteria is not met.
                   """)
        spec.input(
            'relax.convergence_volume',
            valid_type=get_data_class('float'),
            required=False,
            default=get_data_node('float', 0.01),
            help="""
                   The cutoff value for the convergence check on volume. If ``convergence_absolute``
                   is True in AA, otherwise in relative.
                   """)
        spec.input(
            'relax.convergence_positions',
            valid_type=get_data_class('float'),
            required=False,
            default=get_data_node('float', 0.01),
            help="""
                   The cutoff value for the convergence check on positions. If ``convergence_absolute``
                   is True in AA, otherwise in relative difference.
                   """)
        spec.input(
            'relax.convergence_shape_lengths',
            valid_type=get_data_class('float'),
            required=False,
            default=get_data_node('float', 0.1),
            help="""
                   The cutoff value for the convergence check on the lengths of the unit cell
                   vecotrs. If ``convergence_absolute``
                   is True in AA, otherwise in relative difference.
                   """)
        spec.input(
            'relax.convergence_shape_angles',
            valid_type=get_data_class('float'),
            required=False,
            default=get_data_node('float', 0.1),
            help="""
                   The cutoff value for the convergence check on the angles of the unit cell.
                   If ``convergence_absolute`` is True in degrees, otherwise in relative difference.
                   """)
        spec.exit_code(0, 'NO_ERROR', message='the sun is shining')
        spec.exit_code(
            300, 'ERROR_MISSING_REQUIRED_OUTPUT', message='the called workchain does not contain the necessary relaxed output structure')
        spec.exit_code(420, 'ERROR_NO_CALLED_WORKCHAIN', message='no called workchain detected')
        spec.exit_code(500, 'ERROR_UNKNOWN', message='unknown error detected in the relax workchain')
        spec.exit_code(502, 'ERROR_OVERRIDE_PARAMETERS', message='there was an error overriding the parameters')
        spec.outline(
            cls.initialize,
            if_(cls.perform_relaxation)(
                while_(cls.run_next_workchains)(
                    cls.init_next_workchain,
                    cls.run_next_workchain,
                    cls.verify_next_workchain,
                    cls.analyze_convergence,
                ),
                cls.store_relaxed,
            ),
            cls.init_relaxed,
            cls.init_next_workchain,
            cls.run_next_workchain,
            cls.verify_next_workchain,
            cls.results,
            cls.finalize
        )  # yapf: disable

        spec.expose_outputs(cls._next_workchain)
        spec.output('structure_relaxed', valid_type=get_data_class('structure'), required=False)

    def _set_ibrion(self, parameters):
        """Set the algorithm to use for relaxation."""
        if self.inputs.relax.positions.value or self.inputs.relax.shape.value or self.inputs.relax.volume.value:
            parameters.ibrion = self.AlgoEnum.IONIC_RELAXATION_CG
        else:
            parameters.ibrion = self.AlgoEnum.NO_UPDATE

    def _set_ediffg(self, parameters):
        """Set the cutoff to use for relaxation."""
        energy_cutoff = False
        try:
            parameters.ediffg = self.inputs.relax.energy_cutoff.value
            energy_cutoff = True
        except AttributeError:
            pass
        try:
            parameters.ediffg = -abs(self.inputs.relax.force_cutoff.value)
            if energy_cutoff:
                self.report('User supplied both a force and an energy cutoff for the relaxation. Utilizing the force cutoff.')
        except AttributeError:
            pass

    def _set_nsw(self, parameters):
        """Set the number of ionic steps to perform."""
        parameters.nsw = self.inputs.relax.steps.value

    def _set_isif(self, parameters):
        """Set relaxation mode according to the chosen degrees of freedom."""
        parameters.isif = self.ModeEnum.get_from_dof(
            positions=self.inputs.relax.positions.value, shape=self.inputs.relax.shape.value, volume=self.inputs.relax.volume.value)

    def _add_overrides(self, parameters):
        """Add parameters tag overrides, except the ones controlled by other inputs (for provenance)."""
        overrides = AttributeDict({k.lower(): v for k, v in self.inputs.relax.parameters.get_dict().items()})
        parameters.update(overrides)

    def _init_parameters(self):
        """Set parameters parameters based on other inputs."""
        parameters = AttributeDict()
        try:
            input_parameters = self.inputs.parameters.get_dict()
            check_parameters_relax_entries(input_parameters)
            parameters.update(input_parameters)
        except AttributeError:
            pass
        if self.perform_relaxation():
            if 'parameters' in self.inputs.relax:
                # Add override parameters (user force what to use)
                try:
                    self._add_overrides(parameters)
                except ValueError as err:
                    return compose_exit_code(self.exit_code.ERROR_OVERRIDE_PARAMETERS, str(err))
            else:
                # Add plugin controlled flags
                self._set_ibrion(parameters)
                self._set_isif(parameters)
                self._set_nsw(parameters)
                self._set_ediffg(parameters)

        return parameters

    def initialize(self):
        """Initialize."""
        self._init_context()
        self._init_inputs()
        self._init_structure()
        self._init_settings()

    def _init_context(self):
        """Store exposed inputs in the context."""
        self.ctx.exit_code = self.exit_codes.ERROR_UNKNOWN  # pylint: disable=no-member
        self.ctx.is_converged = False
        self.ctx.iteration = 0
        self.ctx.workchains = []
        self.ctx.inputs = AttributeDict()

    def _init_structure(self):
        """Initialize the structure."""
        self.ctx.current_structure = self.inputs.structure

    def _init_settings(self):
        """Initialize the settings."""
        # Make sure we parse the output structure when we want to perform
        # relaxations (override if contrary entry exists).
        if 'settings' in self.inputs:
            settings = AttributeDict(self.inputs.settings.get_dict())
        else:
            settings = AttributeDict({'parser_settings': {}})
        if self.perform_relaxation():
            dict_entry = {'add_structure': True}
            try:
                settings.parser_settings.update(dict_entry)
            except AttributeError:
                settings.parser_settings = dict_entry
        self.ctx.inputs.settings = settings

    def _init_inputs(self):
        """Initialize the inputs."""
        self.ctx.inputs.parameters = self._init_parameters()
        try:
            self._verbose = self.inputs.verbose.value
        except AttributeError:
            pass

    def _set_default_relax_settings(self):
        """Set default settings."""

    def run_next_workchains(self):
        within_max_iterations = bool(self.ctx.iteration < self.inputs.relax.convergence_max_iterations)
        return bool(within_max_iterations and not self.ctx.is_converged)

    def init_relaxed(self):
        """Initialize a calculation based on a relaxed or assumed relaxed structure."""
        if not self.perform_relaxation():
            if self._verbose:
                self.report('skipping structure relaxation and forwarding input/output to the next workchain.')
        else:
            # For the final static run we do not need to parse the output structure, which
            # is at this point enabled.
            settings = AttributeDict(self.ctx.inputs.settings.get_dict())
            settings.parser_settings['add_structure'] = False
            self.ctx.inputs.settings = settings
            self.ctx.inputs.parameters = self.inputs.parameters
            if self._verbose:
                self.report('performing a final calculation using the relaxed structure.')

    def init_next_workchain(self):
        """Initialize the next workchain calculation."""

        if not self.ctx.is_converged:
            self.ctx.iteration += 1

        try:
            self.ctx.inputs
        except AttributeError:
            raise ValueError('no input dictionary was defined in self.ctx.inputs')

        # Set structure
        self.ctx.inputs.structure = self.ctx.current_structure

        # Add exposed inputs
        self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))

        # Make sure we do not have any floating dict (convert to Dict)
        self.ctx.inputs = prepare_process_inputs(self.ctx.inputs)

    def run_next_workchain(self):
        """Run the next workchain."""
        inputs = self.ctx.inputs
        running = self.submit(self._next_workchain, **inputs)

        if not self.ctx.is_converged and self.perform_relaxation():
            self.report('launching {}<{}> iteration #{}'.format(self._next_workchain.__name__, running.pk, self.ctx.iteration))
        else:
            self.report('launching {}<{}> '.format(self._next_workchain.__name__, running.pk))

        return self.to_context(workchains=append_(running))

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

    def analyze_convergence(self):
        """
        Analyze the convergence of the relaxation.

        Compare the input and output structures of the most recent relaxation run. If volume,
        shape and ion positions are all within a given threshold, consider the relaxation converged.
        """
        workchain = self.ctx.workchains[-1]
        # Double check presence of structure
        if 'structure' not in workchain.outputs:
            self.report('The {}<{}> for the relaxation run did not have an '
                        'output structure and most likely failed. However, '
                        'its exit status was zero.'.format(workchain.__class__.__name__, workchain.pk))
            return self.exit_codes.ERROR_NO_RELAXED_STRUCTURE  # pylint: disable=no-member

        self.ctx.previous_structure = self.ctx.current_structure
        self.ctx.current_structure = workchain.outputs.structure

        converged = True
        if self.inputs.relax.convergence_on.value:
            if self._verbose:
                self.report('Checking the convergence of the relaxation.')
            comparison = compare_structures(self.ctx.previous_structure, self.ctx.current_structure)
            delta = comparison.absolute if self.inputs.relax.convergence_absolute.value else comparison.relative
            if self.inputs.relax.positions.value:
                converged &= self.check_positions_convergence(delta)
            if self.inputs.relax.volume.value:
                converged &= self.check_volume_convergence(delta)
            if self.inputs.relax.shape.value:
                converged &= self.check_shape_convergence(delta)

        if not converged:
            self.ctx.current_restart_folder = workchain.outputs.remote_folder
            if self._verbose:
                self.report('Relaxation did not converge, restarting the relaxation.')
        else:
            if self._verbose:
                self.report('Relaxation is converged, finishing with a final static calculation.')
            self.ctx.is_converged = converged

        return self.exit_codes.NO_ERROR  # pylint: disable=no-member

    def check_shape_convergence(self, delta):
        """Check the difference in cell shape before / after the last iteratio against a tolerance."""
        lengths_converged = bool(delta.cell_lengths.max() <= self.inputs.relax.convergence_shape_lengths.value)
        if not lengths_converged:
            self.report('cell lengths changed by max {}, tolerance is {}'.format(delta.cell_lengths.max(),
                                                                                 self.inputs.relax.convergence_shape_lengths.value))

        angles_converged = bool(delta.cell_angles.max() <= self.inputs.relax.convergence_shape_angles.value)
        if not angles_converged:
            self.report('cell angles changed by max {}, tolerance is {}'.format(delta.cell_angles.max(),
                                                                                self.inputs.relax.convergence_shape_angles.value))

        return bool(lengths_converged and angles_converged)

    def check_volume_convergence(self, delta):
        volume_converged = bool(delta.volume <= self.inputs.relax.convergence_volume.value)
        if not volume_converged:
            self.report('cell volume changed by {}, tolerance is {}'.format(delta.volume, self.inputs.relax.convergence_volume.value))
        return volume_converged

    def check_positions_convergence(self, delta):
        positions_converged = bool(delta.pos_lengths.max() <= self.inputs.relax.convergence_positions.value)
        if not positions_converged:
            self.report('max site position change is {}, tolerance is {}'.format(delta.pos_lengths.max(),
                                                                                 self.inputs.relax.convergence_positions.value))
        return positions_converged

    def store_relaxed(self):
        """Store the relaxed structure."""
        workchain = self.ctx.workchains[-1]

        relaxed_structure = workchain.outputs.structure
        if self._verbose:
            self.report("attaching the node {}<{}> as '{}'".format(relaxed_structure.__class__.__name__, relaxed_structure.pk,
                                                                   'structure_relaxed'))
        self.out('structure_relaxed', relaxed_structure)

    def results(self):
        """Attach the remaining output results."""

        workchain = self.ctx.workchains[-1]
        self.out_many(self.exposed_outputs(workchain, self._next_workchain))

    def finalize(self):
        """Finalize the workchain."""

    def perform_relaxation(self):
        """Check if a relaxation is to be performed."""
        return self.inputs.relax.perform.value


def check_parameters_relax_entries(parameters):
    """Check that some relaxation flags are not present in the parameters (no override is allowed)."""

    overrides = AttributeDict({k.lower(): v for k, v in parameters.items()})
    if 'ibrion' in overrides:
        raise ValueError('overriding IBRION not allowed, use inputs to control')
    if 'isif' in overrides:
        raise ValueError('overriding ISIF not allowed, use inputs to control')
    if 'nsw' in overrides:
        raise ValueError('overriding NSW not allowed, use inputs to control')
