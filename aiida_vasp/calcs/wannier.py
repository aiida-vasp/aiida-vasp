# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""Wannier90 - Calculation: Generic"""
from aiida.common.datastructures import CalcInfo, CodeInfo
from aiida.orm import JobCalculation, DataFactory

from aiida_vasp.calcs.base import CalcMeta, Input
from aiida_vasp.calcs.w90win import DictToWin


class WannierBase(object):
    """Base class for Wannier calculations"""

    def write_win(self, inputdict, dst):
        with open(dst, 'w') as win:
            win.write(DictToWin.parse(inputdict[self._win_lname()].get_dict()))

    @staticmethod
    def new_wannier_parameters(*args, **kwargs):
        return DataFactory('parameter')(*args, **kwargs)

    @staticmethod
    def new_wannier_data(*args, **kwargs):
        return DataFactory('vasp.archive')(*args, **kwargs)

    def write_wdat(self, inputdict, path):
        inputdict[self._wdat_lname()].archive.extractall(path)

    @staticmethod
    def _win_lname():
        return 'parameters'

    @staticmethod
    def _wdat_lname():
        return 'data'


class WannierCalculation(JobCalculation, WannierBase):  # pylint: disable=invalid-metaclass
    """Generic Wannier90 calculation, deprecated - to be replaced by aiidateam/aiida-wannier"""
    __metaclass__ = CalcMeta
    imput_file_name = 'wannier90.win'
    output_file_name = 'wannier90.wout'
    parameters = Input(types='parameter')
    data = Input(types='vasp.archive')
    default_parser = 'vasp.wannier'

    def _prepare_for_submission(self, tempfolder, inputdict):
        ''' write input files, set retrieve_list to all files starting with wannier90'''
        # inputfile destinations
        win_dst = tempfolder.get_abs_path('wannier90.win')
        fpath = tempfolder.get_abs_path('')
        # check inputs
        self.verify_inputs(inputdict)
        # write inputs
        self.write_win(inputdict, win_dst)
        self.write_wdat(inputdict, fpath)
        # create and return calcinfo
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.retrieve_list = [['wannier90*', '.', 0]]
        codeinfo = CodeInfo()
        codeinfo.code_uuid = self.get_code().uuid
        codeinfo.code_pk = self.get_code().pk
        codeinfo.cmdline_params = ['wannier90']
        codeinfo.withmpi = False
        calcinfo.codes_info = [codeinfo]
        return calcinfo

    def _init_internal_params(self):
        super(WannierCalculation, self)._init_internal_params()
        self._update_internal_params()

    def verify_inputs(self, inputdict):  # pylint: disable=no-self-use
        '''make sure parameters and data nodes were given'''
        msg = 'input not set: %s'
        if not inputdict.get('parameters'):
            raise ValueError(msg % 'parameters')
        elif not inputdict.get('data'):
            raise ValueError(msg % 'data')
        return True
