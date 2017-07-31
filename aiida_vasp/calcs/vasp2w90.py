from base import Input
from vasp import VaspCalculation
from aiida.orm import DataFactory
from aiida.orm.data.base import List
from aiida_wannier90.io import write_win


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


class Vasp2w90Calculation(VaspCalculation):
    '''General purpose Calculation for using vasp with the vasp2wannier90 interface.'''

    default_parser = 'vasp.vasp2w90'
    wannier_parameters = Input(
        types=['parameter'],
        doc='parameter node: parameters for the wannier interface'
    )
    wannier_projections = Input(
        types=['orbital', List],
        doc='Projections to be defined in the Wannier90 input file.'
    )
    _DEFAULT_PARAMETERS = {'lwannier90': True}

    def write_win(self, inputdict, dst):
        '''Write Wannier90 input file'''
        write_win(
            filename=dst,
            parameters=inputdict.get('wannier90_parameters', {}),
            structure=inputdict.get('structure', None),
            projections=inputdict.get('projections', None)
        )

    def new_wannier_parameters(self, **kwargs):
        return DataFactory('parameter')(**kwargs)

    def write_additional(self, tempfolder, inputdict):
        super(Vasp2w90Calculation, self).write_additional(
            tempfolder, inputdict)
        win = tempfolder.get_abs_path('wannier90.win')
        self.write_win(inputdict, win)
