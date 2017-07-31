__doc__ = '''
This module contains the base class for other vasp parsers.
'''
import re


class BaseParser(object):
    empty_line = re.compile(r'[\r\n]\s*[\r\n]')

    @classmethod
    def line(cls, fobj_or_str, dt=str):
        if isinstance(fobj_or_str, str):
            line = fobj_or_str
        else:
            line = fobj_or_str.readline()
        res = map(dt, line.split())
        if len(res) == 1:
            return res[0]
        else:
            return res

    @classmethod
    def splitlines(cls, fobj_or_str, dt=float):
        if isinstance(fobj_or_str, str):
            lines = fobj_or_str.split('\n')
        else:
            lines = fobj_or_str.readlines()
        return map(lambda l: cls.line(l, dt), lines)


class KeyValueParser(BaseParser):
    '''
    contains regex and functions to find grammar elements
    in vasp input and output files
    '''
    assignment = re.compile(r'(\w+)\s*[=: ]\s*([^;\n]*);?')
    comments = True

    @classmethod
    def get_lines(cls, filename):
        with open(filename) as f:
            lines = f.read().splitlines()
        return lines

    @classmethod
    def retval(cls, *args, **kwargs):
        ret = list(args)
        if cls.comments:
            ret.append(kwargs.get('comment', ''))
        return tuple(ret)

    @classmethod
    def flatten(cls, lst):
        return [i for j in lst for i in j]

    @classmethod
    def find_kv(cls, line):
        return re.findall(cls.assignment, line)

    @classmethod
    def float(cls, string_):
        vals = string_.split()
        value = float(vals.pop(0))
        comment = ' '.join(vals)
        return cls.retval(value, comment=comment)

    @classmethod
    def float_unit(cls, string_):
        vals = string_.split()
        value = float(vals.pop(0))
        unit = vals and vals.pop(0) or ''
        comment = ' '.join(vals)
        return cls.retval(value, unit, comment=comment)

    @classmethod
    def string(cls, string_):
        vals = string_.split()
        value = vals.pop(0)
        comment = ' '.join(vals)
        return cls.retval(value, comment=comment)

    @classmethod
    def bool(cls, string_):
        vals = string_.split()
        bool_str = vals.pop(0)
        if bool_str == 'T':
            value = True
        elif bool_str == 'F':
            value = False
        else:
            raise ValueError('bool string must be "F" or "T"')
        comment = ' '.join(vals)
        return cls.retval(value, comment=comment)

    @classmethod
    def kv_list(cls, filename):
        with open(filename) as potcar:
            kv_list = filter(None, map(cls.find_kv, potcar))
        return kv_list

    @classmethod
    def kv_dict(cls, kv_list):
        kv_dict = dict(cls.flatten(kv_list))
        return kv_dict
