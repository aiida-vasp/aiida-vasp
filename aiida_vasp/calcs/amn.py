# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""Vasp2Wannier90 Calculation: win & mmn -> amn"""
from aiida_vasp.calcs.base import Input
from aiida_vasp.calcs.wannier import WannierBase
from aiida_vasp.calcs.nscf import NscfCalculation


class AmnCalculation(WannierBase, NscfCalculation):
    '''Calculation for using vasp 5 with the vasp2wannier90 interface with
    wannier90 input parameters.'''

    default_parser = 'vasp.amn'
    wannier_settings = Input(
        types=['parameter'],
        doc='parameter node: settings for the ' + 'wannier interface')
    wannier_data = Input(types=['vasp.archive'])

    def _prepare_for_submission(self, tempfolder, inputdict):
        '''remove EIGENVAL and DOSCAR from retrieve_list, so only basc output and
        wannier90 output files are retrieved'''
        calcinfo = super(AmnCalculation, self)._prepare_for_submission(
            tempfolder, inputdict)
        calcinfo.retrieve_list.remove('EIGENVAL')
        calcinfo.retrieve_list.remove('DOSCAR')
        return calcinfo

    def write_additional(self, tempfolder, inputdict):
        super(AmnCalculation, self).write_additional(tempfolder, inputdict)
        windst = tempfolder.get_abs_path('wannier90.win')
        self.write_win(inputdict, windst)
        if inputdict.get('wannier_data'):
            path = tempfolder.get_abs_path('')
            self.write_wdat(inputdict, path)

    def _win_lname(self):
        return 'wannier_settings'

    def _wdat_lname(self):
        return 'wannier_data'
