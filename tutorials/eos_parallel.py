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
from scipy.interpolate import interp1d
from scipy.optimize import minimize

from aiida.common.extendeddicts import AttributeDict
from aiida.engine import WorkChain, append_, calcfunction, while_
from aiida.plugins import DataFactory, WorkflowFactory

from aiida_vasp.utils.workchains import compose_exit_code, prepare_process_inputs


class EosParallelWorkChain(WorkChain):
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
        super(EosParallelWorkChain, cls).define(spec)
        spec.expose_inputs(cls._next_workchain, exclude=['structure'])
        spec.input_namespace(
            'structures', valid_type=DataFactory('structure'), dynamic=True, help='a dictionary of structures to use'
        )
        spec.exit_code(0, 'NO_ERROR', message='the sun is shining')
        spec.exit_code(420, 'ERROR_NO_CALLED_WORKCHAIN', message='no called workchain detected')
        spec.exit_code(500, 'ERROR_UNKNOWN', message='unknown error detected in the eos workchain')

        spec.outline(
            cls.initialize,
            cls.init_and_run_next_workchains,
            cls.verify_next_workchains,
            cls.extract_volume_and_energy,
            cls.finalize
        )  # yapf: disable

        spec.output(
            'eos', valid_type=DataFactory('array'), help='a list containing the cell volumes and total energies'
        )
        spec.output(
            'eos_minimum',
            valid_type=DataFactory('dict'),
            help='a dictionary containing the cell volume at energy minimum'
        )

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

    def init_and_run_next_workchains(self):
        """Initialize and run the workchains."""

        # Check that the context inputs exists
        try:
            self.ctx.inputs
        except AttributeError:
            raise ValueError('No input dictionary was defined in self.ctx.inputs')

        # Take the exposed inputs and add them to the context input. This would
        # be the inputs you supply to this workchain
        self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))

        for structure_key, structure in self.ctx.structures.items():
            print(structure_key, structure)
            # Elevate iteration index
            self.ctx.iteration += 1
            # Take this structure
            self.ctx.inputs.structure = structure
            # Make sure we do not have any floating dict (convert to Dict etc.)
            self.ctx.inputs = prepare_process_inputs(self.ctx.inputs, namespaces=['dynamics'])
            # Submit a VaspWorkChain for this structure
            running = self.submit(self._next_workchain, **self.ctx.inputs)
            self.report(f'launching {self._next_workchain.__name__}<{running.pk}> iteration #{self.ctx.iteration}')
            # Put it into the context and continue
            self.to_context(workchains=append_(running))

    def verify_next_workchains(self):
        """Correct for unexpected behavior."""

        try:
            workchain = self.ctx.workchains[-1]
        except IndexError:
            self.report(f'There is no {self._next_workchain.__name__} in the called workchain list.')
            return self.exit_codes.ERROR_NO_CALLED_WORKCHAIN  # pylint: disable=no-member

        for workchain in self.ctx.workchains:
            # Inherit exit status from last workchain (supposed to be
            # successfull)
            next_workchain_exit_status = workchain.exit_status
            next_workchain_exit_message = workchain.exit_message
            if not next_workchain_exit_status:
                self.ctx.exit_code = self.exit_codes.NO_ERROR  # pylint: disable=no-member
            else:
                self.ctx.exit_code = compose_exit_code(next_workchain_exit_status, next_workchain_exit_message)
                self.report(
                    'The called {}<{}> returned a non-zero exit status. '
                    'The exit status {} is inherited'.format(
                        workchain.__class__.__name__, workchain.pk, self.ctx.exit_code
                    )
                )

        return self.ctx.exit_code

    def extract_volume_and_energy(self):
        """Extract the cell volume and total energy for this structure."""

        for workchain in self.ctx.workchains:
            # Fetch the total energy
            misc = workchain.outputs.misc.get_dict()
            total_energy = misc['total_energies']['energy_extrapolated']

            # Fetch the volume
            volume = workchain.inputs.structure.get_cell_volume()

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
