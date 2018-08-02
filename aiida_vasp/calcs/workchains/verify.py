"""
Base WorkChain for VASP, Error Handling enriched wrapper around VaspCalculation.

Intended to be reused (launched instead of a VaspCalculation) in all other VASP workchains.
Any validation and / or error handling that applies to *every* VASP run,
should be handled on this level, so that every workchain can profit from it.
Anything related to a subset of use cases must be handled in a subclass.
"""
from aiida.work.workchain import WorkChain, while_, append_
from aiida.orm import WorkflowFactory, Code
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
        spec.input('potcar_family', valid_type=get_data_class('str'))
        spec.input('potcar_mapping', valid_type=get_data_class('parameter'))
        spec.input('incar', valid_type=get_data_class('parameter'))
        spec.input('settings', valid_type=get_data_class('parameter'), required=False)
        spec.input('options', valid_type=get_data_class('parameter'))

        spec.outline(
            cls.init_context,
            while_(cls.run_next_workchains)(
                cls.run_next_workchain,
                cls.verify_next_workchain
            ),
            cls.results
        )  # yapf: disable

        spec.output('output_parameters', valid_type=get_data_class('parameter'))
        spec.output('remote_folder', valid_type=get_data_class('remote'))
        spec.output('retrieved', valid_type=get_data_class('folder'))
        spec.output('output_band', valid_type=get_data_class('array.bands'), required=False)
        spec.output('output_structure', valid_type=get_data_class('structure'), required=False)
        spec.output('output_kpoints', valid_type=get_data_class('array.kpoints'), required=False)

    def init_context(self):
        """Initialize context variables that are used during the logical flow of the BaseRestartWorkChain."""
        self.ctx.is_finished = False
        self.ctx.iteration = 0
        self.ctx.max_iterations = self.inputs.max_iterations.value
        self.ctx.inputs = init_input(self.inputs)
        return

    def run_next_workchains(self):
        """
        Return whether a new calculation should be run.

        This is the case as long as the last calculation has not finished successfully and the maximum number of restarts
        has not yet been exceeded.
        """
        return not self.ctx.is_finished and self.ctx.iteration < self.ctx.max_iterations

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
        return self.ctx.workchains[-1].exit_status

    def results(self):
        """Attach the outputs specified in the output specification from the last completed calculation."""
        self.report('workchain completed after {} iterations'.format(self.ctx.iteration))

        workchain = self.ctx.workchains[-1]

        for name, port in self.spec().outputs.iteritems():
            if port.required and name not in workchain.out:
                self.report('the spec specifies the output {} as required but was not an output of {}<{}>'.format(
                    name, self._next_workchain.__name__, workchain.pk))

            if name in workchain.out:
                node = workchain.out[name]
                self.out(name, workchain.out[name])
                if self._verbose:
                    self.report("attaching the node {}<{}> as '{}'".format(node.__class__.__name__, node.pk, name))

        return
