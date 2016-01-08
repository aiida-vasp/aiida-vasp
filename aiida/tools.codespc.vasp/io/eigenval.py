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

    @classmethod
    def parse_eigenval(cls, filename):
        with open(filename) as eig:
            l0 = eig.readline().split()             # read header
            l1 = eig.readline().split()             # "
            l2 = eig.readline().strip()             # "
            l3 = eig.readline().strip()             # "
            name = eig.readline().strip()           # read name line (can be empty)
            p1, nkp, nb = map(int, eig.readline().split()) # read: ? #kp #bands
            data = eig.read()                       # rest is data
        data = re.split(cls.empty_line, data)       # list of data blocks
        data = map(lambda s: s.splitlines(), data)  # list of list of lines
        data = map(lambda s: map(lambda ss: ss.split(), s), data) # 3d list of numbers
        kp = np.zeros((nkp, 4))
        bs = np.zeros((nkp, nb))
        for k, field in enumerate(data):            # iterate over data blocks
            kpbs = filter(None, field)              # throw away empty lines
            kpi = map(float, kpbs.pop(0))           # first line of block is kpoints -> pop
            kp[k] = kpi
            for p in kpbs:                          # rest are band energies
                bs[k, int(p[0])-1] = p[1]           # place energy value in bs[kp, nb] (BandstrucureData format)
        header = {}                                 # build header dict
        header[0] = map(int, l0)
        header[1] = map(float, l1)
        header[2] = float(l2)
        header[3] = l3
        header['name'] = name
        header['some_num'] = p1
        header['n_bands'] = nb
        header['n_kp'] = nkp

        return header, kp, bs
