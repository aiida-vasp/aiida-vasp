# pylint: disable=attribute-defined-outside-init
"""
VerifyWorkChain.

Indented to be used to verify a calculation, perform corrections in inputs files and
restart depending on physical principles etc. E.g. issues that are outside the Calculators awereness,
or not currently checked in it.
"""
from aiida.work.workchain import WorkChain, while_, append_
from aiida.orm import WorkflowFactory, Code
from aiida.orm.data.base import Int, Bool
from aiida_vasp.utils.aiida_utils import get_data_class, init_input


class VerifyWorkChain(WorkChain):
    """Verify the calculations based on basic principles from physics, chemistry and material science."""

    _verbose = True
    _next_workchain = WorkflowFactory('vasp.vasp')

    @classmethod
    def define(cls, spec):
        super(VerifyWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=(get_data_class('structure'), get_data_class('cif')))
        spec.input('kpoints', valid_type=get_data_class('array.kpoints'))
        spec.input('potential_family', valid_type=get_data_class('str'))
        spec.input('potential_mapping', valid_type=get_data_class('parameter'))
        spec.input('incar', valid_type=get_data_class('parameter'))
        spec.input('options', valid_type=get_data_class('parameter'))
        spec.input('settings', valid_type=get_data_class('parameter'), required=False)
        spec.input('restart.max_iterations', valid_type=get_data_class('int'), required=False)
        spec.input('restart.clean_workdir', valid_type=get_data_class('bool'), required=False)
        spec.input(
            'verify.max_iterations',
            valid_type=Int,
            default=Int(1),
            required=False,
            help="""
            the maximum number of iterations VerifyWorkChain will attempt to get the calculation to finish successfully
            """)
        spec.input(
            'verify.clean_workdir',
            valid_type=Bool,
            default=Bool(False),
            required=False,
            help="""
            when set to True, the work directories of all called calculation will be cleaned at the end of VerifyWorkChain execution
            """)
        spec.outline(
            cls.initialize,
            while_(cls.run_next_workchains)(
                cls.run_next_workchain,
                cls.verify_next_workchain
            ),
            cls.results,
            cls.finalize
        )  # yapf: disable

        spec.output('output_parameters', valid_type=get_data_class('parameter'))
        spec.output('remote_folder', valid_type=get_data_class('remote'))
        spec.output('retrieved', valid_type=get_data_class('folder'))
        spec.output('output_structure', valid_type=get_data_class('structure'))

    def initialize(self):
        """Initialize."""
        self._init_context()
        self._init_inputs()

        return

    def _init_context(self):
        """Initialize context variables that are used during the logical flow of the BaseRestartWorkChain."""
        self.ctx.is_finished = False
        self.ctx.iteration = 0
        self.ctx.max_iterations = self.inputs.verify.max_iterations.value

        return

    def _init_inputs(self):
        """Initialize inputs."""
        self.ctx.inputs = init_input(self.inputs, exclude='verify')

        return

    def run_next_workchains(self):
        """
        Return whether a new calculation should be run.

        This is the case as long as the last calculation has not finished successfully and the maximum number of restarts
        has not yet been exceeded.
        """
        return not self.ctx.is_finished and self.ctx.iteration <= self.ctx.max_iterations

    def run_next_workchain(self):
        """Run the lower level VASP workchain."""

        self.ctx.iteration += 1

        try:
            inputs = self.ctx.inputs
        except AttributeError:
            raise ValueError('No input dictionary was defined in self.ctx.inputs')

        running = self.submit(self._next_workchain, **inputs)

        if hasattr(running, 'pid'):
            self.report('launching {}<{}> iteration #{}'.format(self._next_workchain.__name__, running.pid, self.ctx.iteration))
        else:
            # Aiida < 1.0
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

        # Adopt exit status from last child workchain (supposed to be successfull)
        next_workchain = self.ctx.workchains[-1]
        next_workchain_exit_status = next_workchain.exit_status
        if not next_workchain_exit_status:
            self.exit_status = 0
            return
        else:
            self.exit_status = next_workchain_exit_status
            self.report('The child {}<{}> returned a non-zero exit status, {}<{}> '
                        'inherits exit status {}'.format(next_workchain.__class__.__name__, next_workchain.pk, self.__class__.__name__,
                                                         self.pid, next_workchain_exit_status))
        return

    def results(self):
        """Attach the outputs specified in the output specification from the last completed calculation."""

        if not self.exit_status:
            self.report('{}<{}> completed after {} iterations'.format(self.__class__.__name__, self.pid, self.ctx.iteration))

            workchain = self.ctx.workchains[-1]

            for name, port in self.spec().outputs.iteritems():
                if port.required and name not in workchain.out:
                    self.report('the spec specifies the output {} as required '
                                'but was not an output of {}<{}>'.format(name, workchain.__class__.__name__, workchain.pk))

                if name in workchain.out:
                    node = workchain.out[name]
                    self.out(name, workchain.out[name])
                    if self._verbose:
                        self.report("attaching the node {}<{}> as '{}'".format(node.__class__.__name__, node.pk, name))

        return

    def finalize(self):
        """Finalize the workchain."""
        return self.exit_status
