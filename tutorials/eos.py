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
from scipy.optimize import minimize
from scipy.interpolate import interp1d

from aiida.common.extendeddicts import AttributeDict
from aiida.engine import calcfunction, WorkChain, while_, append_
from aiida.plugins import WorkflowFactory, DataFactory
from aiida_vasp.utils.workchains import prepare_process_inputs, compose_exit_code


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
                cls.extract_volume_and_energy
            ),
            cls.finalize
        )  # yapf: disable

        spec.output('eos', valid_type=DataFactory('array'), help='a list containing the cell volumes and total energies')
        spec.output('eos_minimum', valid_type=DataFactory('dict'), help='a dictionary containing the cell volume at energy minimum')

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

        # Define the total energies list
        self.ctx.total_energies = []

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

        # Elavate iteraetion index
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

        self.report('launching {}<{}> iteration #{}'.format(self._next_workchain.__name__, running.pk, self.ctx.iteration))
        return self.to_context(workchains=append_(running))

    def verify_next_workchain(self):
        """Correct for unexpected behavior."""

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

        # Stop further execution of workchains if there are no more structure
        # entries in the structures dictionary
        if not self.ctx.structures:
            self.ctx.is_finished = True

        return self.ctx.exit_code

    def extract_volume_and_energy(self):
        """Extract the cell volume and total energy for this structure."""

        workchain = self.ctx.workchains[-1]

        # Fetch the total energy
        misc = workchain.outputs.misc.get_dict()
        total_energy = misc['total_energies']['energy_no_entropy']

        # Fetch the volume
        volume = self.ctx.inputs.structure.get_cell_volume()

        # Store both in a list
        self.ctx.total_energies.append([volume, total_energy])

    def finalize(self):
        """
        Finalize the workchain.

        Take the total energies container and set is as an output of this workchain.
        """
        # Due to data provenance we cannot return AiiDA data containers that have
        # not been passed through a calcfunction, workfunction or a workchain. Create this now.
        total_energies = store_total_energies(DataFactory('list')(list=self.ctx.total_energies))

        # Let us also try to find a better minimum, just as an example using the power of Python
        energy = locate_minimum(total_energies)

        # And then store the output on the workchain
        self.out('eos', total_energies)
        self.out('eos_minimum', energy)


@calcfunction
def store_total_energies(total_energies):
    """Stores the total energies in ArrayData to keep data provenance."""
    total_energies_list = total_energies.get_list()
    # Let us also sort by volume as we picked the entries in the structures by random
    # above
    total_energies_array = np.array(total_energies_list)
    total_energies_array_sorted = total_energies_array[total_energies_array[:, 0].argsort()]
    array_data = DataFactory('array')()
    array_data.set_array('eos', total_energies_array_sorted)

    return array_data


@calcfunction
def locate_minimum(total_energies):
    """Locate the volume with the lowest energy using interpolation."""
    total_energies_array = total_energies.get_array('eos')
    volumes = total_energies_array[:, 0]
    energies = total_energies_array[:, 1]

    # Establish some initial guess (take lowest energy point from the original dataset)
    min_energy_guess_index = np.argmin(energies)

    # Create the function that can be used to extract interpolated values
    # Using cubic interpolation here, which is not necessarly physically correct,
    # only serves as an example. Please have a look at the theory behind more realiastic
    # models that can be used to fit such data.
    new_energies = interp1d(volumes, energies, kind='cubic')

    # Use minimize from scipy to locate the minimal point that can be ejected by the
    # interpolation routines
    minimized_point = minimize(new_energies, volumes[min_energy_guess_index], tol=1e-3)

    # Create a dictionary to house the results and return
    dict_data = DataFactory('dict')(dict={'volume': minimized_point.x[0], 'energy': minimized_point.fun})

    return dict_data
