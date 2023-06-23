"""
VASP to Wannier90 calculation.

------------------------------
VASP2Wannier90 - Calculation.
"""
from aiida.orm import List
# pylint: disable=abstract-method, unreachable, undefined-variable
# explanation: pylint wrongly complains about (aiida) Node not implementing query
from aiida.plugins import DataFactory

from aiida_vasp.calcs.vasp import VaspCalculation
from aiida_vasp.utils.aiida_utils import get_data_class

#from aiida_wannier90.io import write_win  # pylint: disable=wrong-import-order, import-error


class Vasp2w90Calculation(VaspCalculation):
    """General purpose Calculation for using vasp with the vasp2wannier90 interface."""

    _default_parser = 'vasp.vasp2w90'

    @classmethod
    def define(cls, spec):
        super(Vasp2w90Calculation, cls).define(spec)
        spec.input(
            'wannier_parameters',
            valid_type=get_data_class('core.dict'),
            required=False,
            help='Input parameters for the Wannier90 interface.',
        )
        spec.input(
            'wannier_projections',
            valid_type=(get_data_class('core.orbital'), List),
            required=False,
            help='Projections to be defined in the Wannier90 input.',
        )

    def prepare_for_submission(self, folder):
        """Override the method such that we can add the flag that executes Wannier90 in library mode."""
        raise NotImplementedError('The AiiDA-Wannier90 plugin is not yet v2 ready.')
        # Create a new parameters node that contain the flag that turns on the Wannier90 library mode such
        # that the Wannier90 inputs are created. Execution of Wannier90 itself is to be done after the execution
        # of this calculation with the aiida-wannier90 plugins calculation.
        parameters = self.inputs.parameters.get_dict()
        parameters.update({'lwannier90': True})
        self.inputs.parameters = DataFactory('core.dict')(dict=parameters)
        # Check if any Wannier90 parameters are given
        try:
            _ = self.inputs.wannier_parameters
        except AttributeError:
            self.inputs.wannier_parameters = DataFactory('core.dict')(dict={})

        # Then call the super function.
        return super().prepare_for_submission(folder)

    def write_win(self, dst):
        """Write Wannier90 input."""
        raise NotImplementedError('The AiiDA-Wannier90 plugin is not yet v2 ready.')
        write_win(filename=dst, parameters=self.inputs.wannier_parameters, projections=self.inputs.wannier_projections)

    @staticmethod
    def new_wannier_parameters(**kwargs):
        return DataFactory('core.dict')(**kwargs)

    def write_additional(self, folder, calcinfo):
        super().write_additional(folder, calcinfo)
        win = folder.get_abs_path('wannier90.win')
        self.write_win(win)
