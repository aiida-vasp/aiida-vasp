"""Utility to convert raw input data to .win format"""


class DictToWin(object):
    """Format parameters given in a dictionary into Wannier90 .win format"""

    @classmethod
    def _bool(cls, val):
        return 'T' if val else 'F'

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
        res = ['begin ' + name]
        res += cls._value(val)
        res += ['end ' + name]
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
        return str(val)

    @classmethod
    def _item(cls, key, val):
        if isinstance(val, (list, tuple)):
            return cls._block(key, val)
        return [cls._assign(key, cls._value(val))]

    @classmethod
    def parse(cls, in_dict):
        res = []
        for key, value in in_dict.iteritems():
            res += cls._item(key, value)
        return '\n'.join(res)
