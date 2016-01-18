__doc__ = '''
This module contains tools to read kpoints and bands from EIGENVALUE files.
'''

import re
from parser import BaseParser
import numpy as np


class EigParser(BaseParser):
    '''
    contains regex and functions to find grammar elements
    in EIGENVALUE files
    '''
    def __init__(self, filename):
        res = self.parse_eigenval(filename)
        self.header = res[0]
        self.kp = res[1]
        self.bands = res[2]

    @classmethod
    def parse_eigenval(cls, filename):
        with open(filename) as eig:
            l0 = cls.line(eig, int)        # read header
            l1 = cls.line(eig, float)      # "
            l2 = cls.line(eig, float)      # "
            coord_type = cls.line(eig)     # "
            name = cls.line(eig)           # read name line (can be empty)
            print l0, l1, l2, coord_type, name
            p1, nkp, nb = cls.line(eig, int)  # read: ? #kp #bands
            data = eig.read()               # rest is data
        ni, na, p00, p01 = l0
        data = re.split(cls.empty_line, data)       # list of data blocks
        data = map(lambda s: s.splitlines(), data)  # list of list of lines
        data = map(lambda s: map(lambda ss: ss.split(), s), data)  # 3d list of numbers
        kp = np.zeros((nkp, 4))
        bs = np.zeros((nkp, nb))
        for k, field in enumerate(data):    # iterate over data blocks
            kpbs = filter(None, field)      # throw away empty lines
            kpi = map(float, kpbs.pop(0))   # first line of block is kpoints -> pop
            kp[k] = kpi
            for p in kpbs:                  # rest are band energies
                bs[k, int(p[0])-1] = p[1]   # place energy value in bs[kp, nb] (BandstrucureData format)
        header = {}                         # build header dict
        header[0] = l0
        header[1] = l1
        header[2] = l2
        header['n_ions'] = ni
        header['n_atoms'] = na
        header['p00'] = p00
        header['p01'] = p01
        header['cartesian'] = coord_type.startswith(('c', 'C'))
        header['name'] = name
        header['some_num'] = p1
        header['n_bands'] = nb
        header['n_kp'] = nkp

        return header, kp, bs
