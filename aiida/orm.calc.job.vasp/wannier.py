from base import CalcMeta, Input
from vasp2w90 import dict_to_win
from aiida.common.datastructures import CalcInfo, CodeInfo
from aiida.orm import JobCalculation, DataFactory


class WannierCalculation(JobCalculation):
    __metaclass__ = CalcMeta
    imput_file_name = 'wannier90.win'
    output_file_name = 'wannier90.wout'
    win = Input(types='parameter')
    data = Input(types='vasp.archive')
    default_parser = 'vasp.wannier'

    def _prepare_for_submission(self, tempfolder, inputdict):
        # inputfile destinations
        win_dst = tempfolder.get_abs_path('wannier90.win')
        fpath = tempfolder.get_abs_path('')
        # check inputs
        self.verify_inputs(inputdict)
        # write inputs
        self.write_win(win_dst)
        self.write_data(fpath)
        # create and return calcinfo
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.retrieve_list = [['wannier90*', '.', 0]]
        codeinfo = CodeInfo()
        codeinfo.code_uuid = self.get_code().uuid
        codeinfo.code_pk = self.get_code().pk
        calcinfo.codes_info = [codeinfo]
        return calcinfo

    def _init_internal_params(self):
        super(WannierCalculation, self)._init_internal_params()
        self._update_internal_params()

    def verify_inputs(self, inputdict):
        msg = 'input not set: %s'
        if not inputdict.get('settings'):
            raise ValueError(msg % 'settings')
        elif not inputdict.get('data'):
            raise ValueError(msg % 'data')
        return True

    def new_settings(self, *args, **kwargs):
        return DataFactory('parameter')(*args, **kwargs)

    def write_win(self, dst):
        with open(dst, 'w') as win:
            win.write(
                dict_to_win.parse(self.inpt.settings.get_dict())
            )

    def new_data(self, *args, **kwargs):
        return DataFactory('vasp.archive')(*args, **kwargs)

    def write_data(self, path):
        self.inp.data.archive.extractall(path)
