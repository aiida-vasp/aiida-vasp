from base import CalcMeta, Input
# ~ from vasp2w90 import dict_to_win
from aiida.common.datastructures import CalcInfo, CodeInfo
from aiida.orm import JobCalculation, DataFactory


class dict_to_win(object):

    @classmethod
    def _bool(cls, val):
        return val and 'T' or 'F'

    @classmethod
    def _seq(cls, val):
        res = []
        for i in val:
            if not isinstance(i, (list, tuple)):
                line = cls._value(i)
            else:
                line = ' '.join(map(cls._value, i))
            res.append(' ' + line)
        return res

    @classmethod
    def _block(cls, name, val):
        res = ['begin '+name]
        res += cls._value(val)
        res += ['end '+name]
        return res

    @classmethod
    def _assign(cls, key, val):
        return '{} = {}'.format(key, val)

    @classmethod
    def _value(cls, val):
        if isinstance(val, (str, unicode)):
            return val
        elif isinstance(val, bool):
            return cls._bool(val)
        elif isinstance(val, (list, tuple)):
            return cls._seq(val)
        else:
            return str(val)

    @classmethod
    def _item(cls, key, val):
        if isinstance(val, (list, tuple)):
            return cls._block(key, val)
        else:
            return [cls._assign(key, cls._value(val))]

    @classmethod
    def parse(cls, in_dict):
        res = []
        for k, v in in_dict.iteritems():
            res += cls._item(k, v)
        return '\n'.join(res)


class WannierBase(object):
    def write_win(self, inputdict, dst):
        with open(dst, 'w') as win:
            win.write(
                dict_to_win.parse(inputdict[self._win_lname()].get_dict())
            )

    def new_wannier_settings(self, *args, **kwargs):
        return DataFactory('parameter')(*args, **kwargs)

    def new_wannier_data(self, *args, **kwargs):
        return DataFactory('vasp.archive')(*args, **kwargs)

    def write_wdat(self, inputdict, path):
        inputdict[self._wdat_lname()].archive.extractall(path)

    def _win_lname(self):
        return 'settings'

    def _wdat_lname(self):
        return 'data'


class WannierCalculation(JobCalculation, WannierBase):
    __metaclass__ = CalcMeta
    imput_file_name = 'wannier90.win'
    output_file_name = 'wannier90.wout'
    settings = Input(types='parameter')
    data = Input(types='vasp.archive')
    default_parser = 'vasp.wannier'

    def _prepare_for_submission(self, tempfolder, inputdict):
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

    def verify_inputs(self, inputdict):
        msg = 'input not set: %s'
        if not inputdict.get('settings'):
            raise ValueError(msg % 'settings')
        elif not inputdict.get('data'):
            raise ValueError(msg % 'data')
        return True
