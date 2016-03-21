from base import Input
from wannier import WannierBase
from nscf import NscfCalculation
from aiida.orm import DataFactory


class AmnCalculation(WannierBase, NscfCalculation):
    '''Calculation for using vasp 5.3 with the vasp2wannier90 interface.'''

    default_parser = 'vasp.amn'
    wannier_settings = Input(types=['parameter'],
                             doc='parameter node: settings for the ' +
                             'wannier interface')
    wannier_data = Input(types=['vasp.archive'])

    def _prepare_for_submission(self, tempfolder, inputdict):
        '''changes the retrieve_list to retrieve only files needed
        to continue with nonscf runs'''
        calcinfo = super(AmnCalculation, self)._prepare_for_submission(
            tempfolder, inputdict)
        calcinfo.retrieve_list.remove('EIGENVAL')
        calcinfo.retrieve_list.remove('DOSCAR')
        return calcinfo

    def write_additional(self, tempfolder, inputdict):
        super(AmnCalculation, self).write_additional(
            tempfolder, inputdict)
        windst = tempfolder.get_abs_path('wannier90.win')
        self.write_win(inputdict, windst)
        path = tempfolder.get_abs_path('')
        self.write_wdat(inputdict, path)

    def _win_lname(self):
        return 'wannier_settings'

    def _wdat_lname(self):
        return 'wannier_data'
