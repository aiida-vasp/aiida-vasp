"""
Relax workchain.

----------------
Structure relaxation workchain, which performs the regular duties of relaxing a structure.
It has been designed such that calling workchains should try to use human readable
parameters instead of the code dependent variables.
"""
# pylint: disable=attribute-defined-outside-init
import numpy as np

from aiida.common.extendeddicts import AttributeDict
from aiida.common.exceptions import NotExistent
from aiida.engine import WorkChain, append_, while_, if_
from aiida.plugins import WorkflowFactory

from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.utils.workchains import compare_structures, prepare_process_inputs, compose_exit_code
from aiida_vasp.assistant.parameters import inherit_and_merge_parameters


class RelaxWorkChain(WorkChain):
    """Structure relaxation workchain."""
    _verbose = False
    _next_workchain_string = 'vasp.verify'
    _next_workchain = WorkflowFactory(_next_workchain_string)

    @classmethod
    def define(cls, spec):
        super(RelaxWorkChain, cls).define(spec)
        spec.expose_inputs(cls._next_workchain, exclude=('parameters', 'structure', 'settings'))
        spec.input('structure', valid_type=(get_data_class('structure'), get_data_class('cif')))
        spec.input('parameters', valid_type=get_data_class('dict'))
        spec.input('settings', valid_type=get_data_class('dict'), required=False)
        spec.input('relax.perform',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
            If True, perform relaxation.
            """)
        spec.input('relax.algo',
                   valid_type=get_data_class('str'),
                   default=lambda: get_data_node('str', 'cg'),
                   help="""
            The algorithm to use during relaxation.
            """)
        spec.input('relax.energy_cutoff',
                   valid_type=get_data_class('float'),
                   required=False,
                   help="""
            The cutoff that determines when the relaxation procedure is stopped. In this
            case it stops when the total energy between two ionic steps is less than the
            supplied value.
            """)
        spec.input('relax.force_cutoff',
                   valid_type=get_data_class('float'),
                   required=False,
                   help="""
            The cutoff that determines when the relaxation procedure is stopped. In this
            case it stops when all forces are smaller than than the
            supplied value.
            """)
        spec.input('relax.steps',
                   valid_type=get_data_class('int'),
                   required=False,
                   default=lambda: get_data_node('int', 60),
                   help="""
                   The number of relaxation steps to perform (updates to the atomic positions,
            unit cell size or shape).
                   """)
        spec.input('relax.positions',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', True),
                   help="""
                   If True, perform relaxation of the atomic positions.
                   """)
        spec.input('relax.shape',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
                   If True, perform relaxation of the unit cell shape.
                   """)
        spec.input('relax.volume',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
                   If True, perform relaxation of the unit cell volume..
                   """)
        spec.input('relax.convergence_on',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
                   If True, test convergence based on selected criterias set.
                   """)
        spec.input('relax.convergence_absolute',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
                   If True, test convergence based on absolute differences.
                   """)
        spec.input('relax.convergence_max_iterations',
                   valid_type=get_data_class('int'),
                   required=False,
                   default=lambda: get_data_node('int', 5),
                   help="""
                   The number of iterations to perform if the convergence criteria is not met.
                   """)
        spec.input('relax.convergence_volume',
                   valid_type=get_data_class('float'),
                   required=False,
                   default=lambda: get_data_node('float', 0.01),
                   help="""
                   The cutoff value for the convergence check on volume. If ``convergence_absolute``
                   is True in AA, otherwise in relative.
                   """)
        spec.input('relax.convergence_positions',
                   valid_type=get_data_class('float'),
                   required=False,
                   default=lambda: get_data_node('float', 0.01),
                   help="""
                   The cutoff value for the convergence check on positions. If ``convergence_absolute``
                   is True in AA, otherwise in relative difference.
                   """)
        spec.input('relax.convergence_shape_lengths',
                   valid_type=get_data_class('float'),
                   required=False,
                   default=lambda: get_data_node('float', 0.1),
                   help="""
                   The cutoff value for the convergence check on the lengths of the unit cell
                   vecotrs. If ``convergence_absolute``
                   is True in AA, otherwise in relative difference.
                   """)
        spec.input('relax.convergence_shape_angles',
                   valid_type=get_data_class('float'),
                   required=False,
                   default=lambda: get_data_node('float', 0.1),
                   help="""
                   The cutoff value for the convergence check on the angles of the unit cell.
                   If ``convergence_absolute`` is True in degrees, otherwise in relative difference.
                   """)
        spec.exit_code(0, 'NO_ERROR', message='the sun is shining')
        spec.exit_code(300,
                       'ERROR_MISSING_REQUIRED_OUTPUT',
                       message='the called workchain does not contain the necessary relaxed output structure')
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
        spec.output('relax.structure', valid_type=get_data_class('structure'), required=False)

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
        self.ctx.relax = False
        self.ctx.iteration = 0
        self.ctx.workchains = []
        self.ctx.inputs = AttributeDict()

    def _init_inputs(self):
        """Initialize the inputs."""
        self.ctx.inputs.parameters = self._init_parameters()
        try:
            self._verbose = self.inputs.verbose.value
            self.ctx.inputs.verbose = self.inputs.verbose
        except AttributeError:
            pass

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

    def _init_parameters(self):
        """Collect input to the workchain in the relax namespace and put that into the parameters."""

        # At some point we will replace this with possibly input checking using the PortNamespace on
        # a dict parameter type. As such we remove the workchain input parameters as node entities. Much of
        # the following is just a workaround until that is in place in AiiDA core.
        parameters = inherit_and_merge_parameters(self.inputs)

        # Need to save the parameter that controls if relaxation should be performed
        self.ctx.relax = parameters.relax.perform
        if not parameters.relax.perform:
            # Make sure we do not expose the relax namespace in the input parameters (
            # basically setting no code tags related to relaxation unless user overrides)
            del parameters.relax

        return parameters

    def _set_default_relax_settings(self):
        """Set default settings."""

    def run_next_workchains(self):
        within_max_iterations = bool(self.ctx.iteration < self.inputs.relax.convergence_max_iterations)
        return bool(within_max_iterations and not self.ctx.is_converged)

    def init_relaxed(self):
        """Initialize a calculation based on a relaxed or assumed relaxed structure."""
        if not self.perform_relaxation():
            if self._verbose:
                self.report('skipping structure relaxation and forwarding to the next workchain.')
        else:
            # For the final static run we do not need to parse the output structure, which
            # is at this point enabled.
            self.ctx.inputs.settings.parser_settings.add_structure = False
            try:
                self.ctx.inputs.settings.parser_settings.add_structure = False
            except AttributeError:
                pass
            # Remove relaxation parameters if they exist
            try:
                del self.ctx.inputs.parameters.relax
            except AttributeError:
                pass
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

        # Make sure we do not have any floating dict (convert to Dict etc.)
        self.ctx.inputs_ready = prepare_process_inputs(self.ctx.inputs, namespaces=['verify', 'dynamics'])

    def run_next_workchain(self):
        """Run the next workchain."""
        inputs = self.ctx.inputs_ready
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
            # Make sure at the very minimum we attach the misc node (if it exists) that contains notifications and other
            # quantities that can be salvaged
            try:
                self.out('misc', workchain.outputs['misc'])
            except NotExistent:
                pass

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
        if self.ctx.inputs.parameters.relax.convergence_on:
            if self._verbose:
                self.report('Checking the convergence of the relaxation.')
            comparison = compare_structures(self.ctx.previous_structure, self.ctx.current_structure)
            delta = comparison.absolute if self.ctx.inputs.parameters.relax.convergence_absolute else comparison.relative
            if self.ctx.inputs.parameters.relax.positions:
                converged &= self.check_positions_convergence(delta)
            if self.ctx.inputs.parameters.relax.volume:
                converged &= self.check_volume_convergence(delta)
            if self.ctx.inputs.parameters.relax.shape:
                converged &= self.check_shape_convergence(delta)

            if not converged:
                self.ctx.current_restart_folder = workchain.outputs.remote_folder
                if self._verbose:
                    self.report('Relaxation did not converge, restarting the relaxation.')
            else:
                if self._verbose:
                    self.report('Relaxation is considered converged.')

        if converged:
            self.ctx.is_converged = converged

        return self.exit_codes.NO_ERROR  # pylint: disable=no-member

    def check_shape_convergence(self, delta):
        """Check the difference in cell shape before / after the last iteratio against a tolerance."""
        lengths_converged = bool(delta.cell_lengths.max() <= self.ctx.inputs.parameters.relax.convergence_shape_lengths)
        if not lengths_converged:
            self.report('cell lengths changed by max {}, tolerance is {}'.format(
                delta.cell_lengths.max(), self.ctx.inputs.parameters.relax.convergence_shape_lengths))

        angles_converged = bool(delta.cell_angles.max() <= self.ctx.inputs.parameters.relax.convergence_shape_angles)
        if not angles_converged:
            self.report('cell angles changed by max {}, tolerance is {}'.format(delta.cell_angles.max(),
                                                                                self.ctx.inputs.parameters.relax.convergence_shape_angles))

        return bool(lengths_converged and angles_converged)

    def check_volume_convergence(self, delta):
        """Check the convergence of the volume, given a cutoff."""
        volume_converged = bool(delta.volume <= self.ctx.inputs.parameters.relax.convergence_volume)
        if not volume_converged:
            self.report('cell volume changed by {}, tolerance is {}'.format(delta.volume,
                                                                            self.ctx.inputs.parameters.relax.convergence_volume))
        return volume_converged

    def check_positions_convergence(self, delta):
        """Check the convergence of the atomic positions, given a cutoff."""
        try:
            positions_converged = bool(np.nanmax(delta.pos_lengths) <= self.ctx.inputs.parameters.relax.convergence_positions)
        except RuntimeWarning:
            # Here we encountered the case of having one atom centered at the origin, so
            # we do not know if it is converged, so settings it to False
            self.report('there is NaN entries in the relative comparison for '
                        'the positions during relaxation, assuming position is not converged')
            positions_converged = False

        if not positions_converged:
            try:
                self.report('max site position change is {}, tolerance is {}'.format(np.nanmax(
                    delta.pos_lengths), self.ctx.inputs.parameters.relax.convergence_positions))
            except RuntimeWarning:
                pass

        return positions_converged

    def store_relaxed(self):
        """Store the relaxed structure."""
        workchain = self.ctx.workchains[-1]

        relaxed_structure = workchain.outputs.structure
        if self._verbose:
            self.report("attaching the node {}<{}> as '{}'".format(relaxed_structure.__class__.__name__, relaxed_structure.pk,
                                                                   'relax.structure'))
        self.out('relax.structure', relaxed_structure)

    def results(self):
        """Attach the remaining output results."""

        workchain = self.ctx.workchains[-1]
        self.out_many(self.exposed_outputs(workchain, self._next_workchain))

    def finalize(self):
        """Finalize the workchain."""

    def perform_relaxation(self):
        """Check if a relaxation is to be performed."""
        return self.ctx.relax
