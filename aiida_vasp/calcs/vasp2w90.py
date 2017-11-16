# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""VASP2Wannier90 - Calculation"""
from aiida.orm import DataFactory
from aiida.orm.data.base import List
from aiida_wannier90.io import write_win

from aiida_vasp.calcs.base import Input
from aiida_vasp.calcs.vasp import VaspCalculation


class Vasp2w90Calculation(VaspCalculation):
    """General purpose Calculation for using vasp with the vasp2wannier90 interface."""

    default_parser = 'vasp.vasp2w90'
    wannier_parameters = Input(types=['parameter'], doc='parameter node: parameters for the wannier interface')
    wannier_projections = Input(types=['orbital', List], doc='Projections to be defined in the Wannier90 input file.')
    _DEFAULT_PARAMETERS = {'lwannier90': True}

    @staticmethod
    def write_win(inputdict, dst):
        """Write Wannier90 input file"""
        write_win(filename=dst, parameters=inputdict.get('wannier_parameters', {}), projections=inputdict.get('wannier_projections', None))

    @staticmethod
    def new_wannier_parameters(**kwargs):
        return DataFactory('parameter')(**kwargs)

    def write_additional(self, tempfolder, inputdict):
        super(Vasp2w90Calculation, self).write_additional(tempfolder, inputdict)
        win = tempfolder.get_abs_path('wannier90.win')
        self.write_win(inputdict, win)
