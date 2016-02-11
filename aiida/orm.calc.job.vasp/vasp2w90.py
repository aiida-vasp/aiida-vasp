from base import Input
from vasp5 import Vasp5Calculation
from aiida.orm import DataFactory


class dict_to_win(object):

    @classmethod
    def _bool(cls, val):
        return val and 'T' or 'F'

    @classmethod
    def _seq(cls, val):
        return map(' '.join, map(cls._value, val))

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


class Vasp2W90Calculation(Vasp5Calculation):
    '''Calculation for using vasp 5.3 with the vasp2wannier90 interface.'''

    wannier_settings = Input(types=['parameter'],
                             doc='parameter node: settings for the ' +
                             'wannier interface')

    def write_wannier(self, inputdict, dst):
        '''convert dict to wannier input and write to file'''
        if 'wannier_settings' in inputdict:
            with open(dst, 'w') as win:
                win.write(dict_to_win.parse(self.inp.wannier_settings.get_dict()))

    def new_wannier(self, **kwargs):
        return DataFactory('parameter')(**kwargs)

    def write_additional(self, tempfolder, inputdict):
        super(Vasp2W90Calculation, self).write_additional(tempfolder, inputdict)
        win = tempfolder.get_abs_path('wannier90.win')
        self.write_wannier(inputdict, win)
