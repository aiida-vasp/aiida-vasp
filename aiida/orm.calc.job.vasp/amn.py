from base import Input
from nscf import NscfCalculation
from aiida.orm import DataFactory


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


class AmnCalculation(NscfCalculation):
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

    def write_wannier(self, inputdict, dst):
        '''convert dict to wannier input and write to file'''
        if 'wannier_settings' in inputdict:
            with open(dst, 'w') as win:
                win.write(
                    dict_to_win.parse(self.inp.wannier_settings.get_dict()))

    def write_wdat(self, inputdict, dst):
        if 'wannier_data' in inputdict:
            self.inp.wannier_data.archive.extractall(path=dst)

    def new_wannier_settings(self, **kwargs):
        return DataFactory('parameter')(**kwargs)

    def new_wannier_data(self, **kwargs):
        return DataFactory('vasp.archive')(**kwargs)

    def write_additional(self, tempfolder, inputdict):
        super(AmnCalculation, self).write_additional(
            tempfolder, inputdict)
        win = tempfolder.get_abs_path('wannier90.win')
        self.write_wannier(inputdict, win)
        path = tempfolder.get_abs_path('')
        self.write_wdat(inputdict, path)
