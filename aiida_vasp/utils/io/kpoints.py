from parser import BaseParser
import numpy as np


class KpParser(BaseParser):

    def __init__(self, fname):
        self.header, res = self.parse_kp_file(fname)
        self.kpoints = res[:, :3]
        if res.shape[1] == 4:
            self.weights = res[:, 3]
        else:
            self.weights = None
        self.cartesian = self.header['cartesian']

    @classmethod
    def parse_kp(cls, fobj_or_str):
        if isinstance(fobj_or_str, str):
            from StringIO import StringIO
            fobj_or_str = StringIO(fobj_or_str)
        header = {}
        header['name'] = fobj_or_str.readline()
        header['nkp'] = cls.line(fobj_or_str, dt=int)
        header['cartesian'] = not fobj_or_str.readline().startswith(('r', 'R'))
        lines = np.array(cls.splitlines(fobj_or_str))
        return header, lines

    @classmethod
    def parse_kp_file(cls, fname):
        with open(fname) as kp:
            return cls.parse_kp(kp)
