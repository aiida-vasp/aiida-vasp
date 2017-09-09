# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""VASP2Wannier90 - Calculation"""
from aiida.orm import DataFactory

from aiida_vasp.calcs.base import Input
from aiida_vasp.calcs.vasp5 import Vasp5Calculation


class Vasp2w90Calculation(Vasp5Calculation):
    '''General purpose Calculation for using vasp with the vasp2wannier90 interface.
    Same data storage space concerns apply as with :py:class:`Vasp5Calculation`.'''

    default_parser = 'vasp.vasp2w90'
    wannier_parameters = Input(
        types=['parameter'],
        doc='parameter node: parameters for the ' + 'wannier interface')

    def write_win(self, inputdict, dst):
        '''convert dict to wannier input and write to file'''
        if 'wannier_parameters' in inputdict:
            super(Vasp2w90Calculation, self).write_win(inputdict, dst)

    def new_wannier_parameters(self, **kwargs):
        return DataFactory('parameter')(**kwargs)

    def new_wannier_data(self, **kwargs):
        return DataFactory('vasp.archive')(**kwargs)

    def write_additional(self, tempfolder, inputdict):
        super(Vasp2w90Calculation, self).write_additional(
            tempfolder, inputdict)
        win = tempfolder.get_abs_path('wannier90.win')
        self.write_win(inputdict, win)


class Vasp2W90Calculation(Vasp2w90Calculation):
    '''Calculation for using vasp 5.3 with the vasp2wannier90 interface.'''
    pass
