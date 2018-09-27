# pylint: disable=attribute-defined-outside-init
"""
VerifyWorkChain.

Indented to be used to verify a calculation, perform corrections in inputs files and
restart depending on physical principles etc. E.g. issues that are outside the Calculators awereness,
or not currently checked in it.
"""
from aiida.plugins.entry_point import load_entry_point_from_string
from aiida.common.extendeddicts import AttributeDict
from aiida.work.workchain import WorkChain, while_, append_
from aiida.orm import WorkflowFactory
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.utils.workchains import prepare_process_inputs


class VerifyWorkChain(WorkChain):
    """Verify the calculations based on basic principles from physics, chemistry and material science."""

    _verbose = False
    _next_workchain_string = 'vasp.vasp'
    _next_workchain = WorkflowFactory(_next_workchain_string)

    @classmethod
    def define(cls, spec):
        super(VerifyWorkChain, cls).define(spec)
        spec.expose_inputs(
            load_entry_point_from_string('aiida.workflows:' + cls._next_workchain_string), exclude=('structure', 'kpoints', 'parameters'))
        spec.input('structure', valid_type=(get_data_class('structure'), get_data_class('cif')))
        spec.input('kpoints', valid_type=get_data_class('array.kpoints'))
        spec.input('parameters', valid_type=get_data_class('parameter'))
        spec.input(
            'max_iterations',
            valid_type=get_data_class('int'),
            required=False,
            default=get_data_node('int', 5),
            help="""
                   The maximum number of iterations to perform.
                   """)
        spec.input(
            'clean_workdir',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', True),
            help="""
                   If True, clean the work dir upon the completion of a successfull calculation.
                   """)

        spec.outline(
            cls.initialize,
            while_(cls.run_next_workchains)(
                cls.init_next_workchain,
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
        spec.output('output_kpoints', valid_type=get_data_class('array.kpoints'), required=False)
        spec.output('output_trajectory', valid_type=get_data_class('array.trajectory'), required=False)
        spec.output('output_chgcar', valid_type=get_data_class('vasp.chargedensity'), required=False)
        spec.output('output_wavecar', valid_type=get_data_class('vasp.wavefun'), required=False)
        spec.output('output_bands', valid_type=get_data_class('array.bands'), required=False)
        spec.output('output_dos', valid_type=get_data_class('array'), required=False)
        spec.output('output_occupations', valid_type=get_data_class('array'), required=False)
        spec.output('output_energies', valid_type=get_data_class('array'), required=False)
        spec.output('output_projectors', valid_type=get_data_class('array'), required=False)
        spec.output('output_dielectrics', valid_type=get_data_class('array'), required=False)
        spec.output('output_born_charges', valid_type=get_data_class('array'), required=False)
        spec.output('output_hessian', valid_type=get_data_class('array'), required=False)
        spec.output('output_dynmat', valid_type=get_data_class('array'), required=False)
        spec.output('output_final_forces', valid_type=get_data_class('array'), required=False)
        spec.output('output_final_stress', valid_type=get_data_class('array'), required=False)

    def initialize(self):
        """Initialize."""
        self._init_context()
        self._init_inputs()

        return

    def _init_context(self):
        """Initialize context variables that are used during the logical flow of the BaseRestartWorkChain."""
        self.ctx.is_finished = False
        self.ctx.iteration = 0
        self.ctx.inputs = AttributeDict()

        return

    def _init_inputs(self):
        """Initialize inputs."""
        # Set structure, kpoints and parameters (these will be modified in the future)
        self.ctx.inputs.structure = self.inputs.structure
        self.ctx.inputs.kpoints = self.inputs.kpoints
        self.ctx.inputs.parameters = self.inputs.parameters
        try:
            self._verbose = self.inputs.verbose.value
        except AttributeError:
            pass

    def run_next_workchains(self):
        """
        Return whether a new calculation should be run.

        This is the case as long as the last calculation has not finished successfully and the maximum number of restarts
        has not yet been exceeded.
        """
        return not self.ctx.is_finished and self.ctx.iteration <= self.inputs.max_iterations.value

    def init_next_workchain(self):
        """Initialize the next workchain."""

        self.ctx.iteration += 1

        try:
            self.ctx.inputs
        except AttributeError:
            raise ValueError('No input dictionary was defined in self.ctx.inputs')

        # Add exposed inputs
        self.ctx.inputs.update(self.exposed_inputs(load_entry_point_from_string('aiida.workflows:' + self._next_workchain_string)))

        # Make sure we do not have any floating dict (convert to ParameterData)
        self.ctx.inputs = prepare_process_inputs(self.ctx.inputs)

    def run_next_workchain(self):
        """Run the next workchain."""
        inputs = self.ctx.inputs
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

        # Adopt exit status from last child workchain (supposed to be
        # successfull)
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
