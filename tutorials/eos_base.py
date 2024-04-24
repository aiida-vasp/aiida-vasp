"""
The eos workchain

----------------

The workchain will accept a dictionary of structures and extract the
total energies for each structure.
The data is saved and the energy minimum is calculated and stored.
"""
# pylint: disable=attribute-defined-outside-init
import random

import numpy as np
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import WorkChain, append_, calcfunction, while_
from aiida.plugins import DataFactory, WorkflowFactory
from aiida_vasp.utils.workchains import compose_exit_code, prepare_process_inputs


class EosWorkChain(WorkChain):
    """
    The eos workchain

    ----------------

    The workchain will accept a dictionary of structures and extract the
    total energies for each structure.
    The data is saved and the energy minimum is calculated and stored.
    """

    _verbose = False
    _next_workchain_string = 'vasp.vasp'
    _next_workchain = WorkflowFactory(_next_workchain_string)

    @classmethod
    def define(cls, spec):
        super(EosWorkChain, cls).define(spec)
        spec.expose_inputs(cls._next_workchain, exclude=['structure'])
        spec.input_namespace('structures', valid_type=DataFactory('structure'), dynamic=True, help='a dictionary of structures to use')
        spec.exit_code(0, 'NO_ERROR', message='the sun is shining')
        spec.exit_code(420, 'ERROR_NO_CALLED_WORKCHAIN', message='no called workchain detected')
        spec.exit_code(500, 'ERROR_UNKNOWN', message='unknown error detected in the eos workchain')

        spec.outline(
            cls.initialize,
            while_(cls.run_next_workchains)(
                cls.init_next_workchain,
                cls.run_next_workchain,
                cls.verify_next_workchain,
                cls.extract_quantities
            ),
            cls.finalize
        )  # yapf: disable

        spec.output('quantities', valid_type=DataFactory('array'), help='a container for various quantities')
        spec.expose_outputs(cls._next_workchain)

    def initialize(self):
        """Initialize the eos workchain."""
        self._init_context()
        self._init_inputs()

    def _init_context(self):
        """Initialize context variables that are used during the."""
        # Set the exit code to error in case we forget to set it to NO_ERROR
        self.ctx.exit_code = self.exit_codes.ERROR_UNKNOWN  # pylint: disable=no-member

        # Copy structures to context, since we will empty it as we go.
        # Since structures is an input to a workchain we cannot modify it and need to copy.
        self.ctx.structures = dict(self.inputs.structures)

        # Continue to submit workchains until this is True
        self.ctx.is_finished = False

        # Define an interation index
        self.ctx.iteration = 0

        # Define the context inputs
        self.ctx.inputs = AttributeDict()

        # Define container to store quantities that is extracted in each step
        self.ctx.quantities_container = []

    def _init_inputs(self):
        """Initialize inputs."""

        # Set verbose flag (more or less details on the workchain report)
        try:
            self._verbose = self.inputs.verbose.value
        except AttributeError:
            pass

    def run_next_workchains(self):
        """
        Return whether a new workchain should be run.

        This is the case as long as the last workchain has not finished successfully.
        """
        return not self.ctx.is_finished

    def init_next_workchain(self):
        """Initialize the next workchain."""

        # Elevate iteration index
        self.ctx.iteration += 1

        # Check that the context inputs exists
        try:
            self.ctx.inputs
        except AttributeError:
            raise ValueError('No input dictionary was defined in self.ctx.inputs')

        # Take the exposed inputs and add them to the context input. This would
        # be the inputs you supply to this workchain
        self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))

        # We did not expose the structure as we would like to set this
        # from the supplied structures input. Just choose any item in the
        # structures dictionary and asign that to the next run.
        item = random.choice(list(self.ctx.structures.keys()))
        self.ctx.inputs.structure = self.ctx.structures.pop(item)

        # Make sure we do not have any floating dict (convert to Dict etc.)
        self.ctx.inputs = prepare_process_inputs(self.ctx.inputs)

    def run_next_workchain(self):
        """
        Run the next workchain

        It is either submitted to the daemon or run, depending on how you
        run this workchain. The execution method is inherited.
        """
        inputs = self.ctx.inputs
        running = self.submit(self._next_workchain, **inputs)
        self.report(f'launching {self._next_workchain.__name__}<{running.pk}> iteration #{self.ctx.iteration}')
        self.to_context(workchains=append_(running))

    def verify_next_workchain(self):
        """Correct for unexpected behavior."""

        try:
            workchain = self.ctx.workchains[-1]
        except IndexError:
            self.report(f'There is no {self._next_workchain.__name__} in the called workchain list.')
            return self.exit_codes.ERROR_NO_CALLED_WORKCHAIN  # pylint: disable=no-member

        # Inherit exit status from last workchain (supposed to be
        # successfull)
        next_workchain_exit_status = workchain.exit_status
        next_workchain_exit_message = workchain.exit_message
        if not next_workchain_exit_status:
            self.ctx.exit_code = self.exit_codes.NO_ERROR  # pylint: disable=no-member
        else:
            self.ctx.exit_code = compose_exit_code(next_workchain_exit_status, next_workchain_exit_message)
            self.report(
                'The called {}<{}> returned a non-zero exit status. ' 'The exit status {} is inherited'.format(
                    workchain.__class__.__name__, workchain.pk, self.ctx.exit_code
                )
            )

        # Stop further execution of workchains if there are no more structure
        # entries in the structures dictionary
        if not self.ctx.structures:
            self.ctx.is_finished = True

        return self.ctx.exit_code

    def extract_quantities(self):
        """Extract the quantities you want and store them in the container."""

        # What happens in this routine needs modifications by you.
        # Typically you get the output of the called workchain by doing
        # workchain = self.ctx.workchains[-1]
        # some_output_quantity = workchain.outputs.some_output_quantity

        # An example which stores nonsense.
        self.ctx.quantities_container.append([self.ctx.iteration, f'Some quantity for iteration: {self.ctx.iteration}'])

    def finalize(self):
        """
        Finalize the workchain.

        Take the quantity container and set it as an output of this workchain.
        """
        # Due to data provenance we cannot return AiiDA data containers that have
        # not been passed through a calcfunction, workfunction or a workchain. Create this now.
        quantities_container = store_quantities(DataFactory('list')(list=self.ctx.quantities_container))

        # And then store the output on the workchain
        self.out('quantities', quantities_container)


@calcfunction
def store_quantities(quantities_container):
    """Stores the quantities to keep data provenance."""
    quantities_container_list = quantities_container.get_list()
    quantities_container_array = DataFactory('array')()
    quantities_container_array.set_array('quantities', np.array(quantities_container_list))
    return quantities_container_array
