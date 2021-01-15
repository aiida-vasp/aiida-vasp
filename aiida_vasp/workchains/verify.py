"""
Verify workchain.

-----------------
Indented to be used to verify a calculation, perform corrections in inputs files and
restart depending on physical principles etc. E.g. issues that are outside the Calculators awereness,
or not currently checked in it. This workchain does currently nothing.
"""

# pylint: disable=attribute-defined-outside-init
from aiida.common.extendeddicts import AttributeDict
from aiida.common.exceptions import NotExistent
from aiida.engine import WorkChain, while_, append_
from aiida.plugins import WorkflowFactory

from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.utils.workchains import prepare_process_inputs, compose_exit_code


class VerifyWorkChain(WorkChain):
    """Verify the calculations based on basic principles from physics, chemistry and material science."""

    _verbose = False
    _next_workchain_string = 'vasp.vasp'
    _next_workchain = WorkflowFactory(_next_workchain_string)

    @classmethod
    def define(cls, spec):
        super(VerifyWorkChain, cls).define(spec)
        spec.expose_inputs(cls._next_workchain)
        spec.input('verify.max_iterations',
                   valid_type=get_data_class('int'),
                   required=False,
                   default=lambda: get_data_node('int', 1),
                   help="""
                   The maximum number of iterations to perform.
                   """)
        spec.exit_code(0, 'NO_ERROR', message='the sun is shining')
        spec.exit_code(420, 'ERROR_NO_CALLED_WORKCHAIN', message='no called workchain detected')
        spec.exit_code(500, 'ERROR_UNKNOWN', message='unknown error detected in the verify workchain')
        spec.outline(
            cls.initialize,
            while_(cls.run_next_workchains)(
                cls.init_next_workchain,
                cls.run_next_workchain,
                cls.verify_next_workchain
            ),
            cls.finalize
        )  # yapf: disable

        spec.expose_outputs(cls._next_workchain)

    def initialize(self):
        """Initialize."""
        self._init_context()
        self._init_inputs()

    def _init_context(self):
        """Initialize context variables that are used during the logical flow."""
        self.ctx.exit_code = self.exit_codes.ERROR_UNKNOWN  # pylint: disable=no-member
        self.ctx.is_finished = False
        self.ctx.iteration = 0
        self.ctx.inputs = AttributeDict()

    def _init_inputs(self):
        """Initialize inputs."""
        try:
            self._verbose = self.inputs.verbose.value
            self.ctx.inputs.verbose = self.inputs.verbose
        except AttributeError:
            pass

    def run_next_workchains(self):
        """
        Return whether a new calculation should be run.

        This is the case as long as the last calculation has not finished successfully and the maximum number of restarts
        has not yet been exceeded.
        """
        return not self.ctx.is_finished and self.ctx.iteration <= self.inputs.verify.max_iterations.value

    def init_next_workchain(self):
        """Initialize the next workchain."""

        self.ctx.iteration += 1

        try:
            self.ctx.inputs
        except AttributeError:
            raise ValueError('No input dictionary was defined in self.ctx.inputs')

        # Add exposed inputs
        self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))

        # Make sure we do not have any floating dict (convert to Dict)
        self.ctx.inputs = prepare_process_inputs(self.ctx.inputs, namespaces=['dynamics'])

    def run_next_workchain(self):
        """Run the next workchain."""
        inputs = self.ctx.inputs
        running = self.submit(self._next_workchain, **inputs)

        self.report('launching {}<{}> iteration #{}'.format(self._next_workchain.__name__, running.pk, self.ctx.iteration))
        return self.to_context(workchains=append_(running))

    def verify_next_workchain(self):
        """
        Correct for unexpected physics/chemistry/material science behavior.

        Here we should correct all things that voids what we expect from
        physics/chemistry/material science. I.e. things that cannot be corrected for at the
        calculation level (simple restarts etc.).
        """

        # Currently only set to finished on first go.
        self.ctx.is_finished = True

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

    def finalize(self):
        """Finalize the workchain."""
        workchain = self.ctx.workchains[-1]
        self.out_many(self.exposed_outputs(workchain, self._next_workchain))
