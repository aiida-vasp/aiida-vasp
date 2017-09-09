from base import Input
from vasp import VaspCalculation
from aiida.orm import DataFactory
from aiida.orm.data.base import List
from aiida_wannier90.io import write_win


class Vasp2w90Calculation(VaspCalculation):
    '''General purpose Calculation for using vasp with the vasp2wannier90 interface.'''

    default_parser = 'vasp.vasp2w90'
    wannier_parameters = Input(
        types=['parameter'],
        doc='parameter node: parameters for the wannier interface')
    wannier_projections = Input(
        types=['orbital', List],
        doc='Projections to be defined in the Wannier90 input file.')
    _DEFAULT_PARAMETERS = {
        'lwannier90': True,
        'ncore': 1,
        'isym': -1,
        'icharg': 11,
        'lpead': False,
        'lwave': False
    }

    def write_win(self, inputdict, dst):
        '''Write Wannier90 input file'''
        write_win(
            filename=dst,
            parameters=inputdict.get('wannier90_parameters', {}),
            projections=inputdict.get('projections', None))

    def new_wannier_parameters(self, **kwargs):
        return DataFactory('parameter')(**kwargs)

    def write_additional(self, tempfolder, inputdict):
        super(Vasp2w90Calculation, self).write_additional(
            tempfolder, inputdict)
        win = tempfolder.get_abs_path('wannier90.win')
        self.write_win(inputdict, win)
