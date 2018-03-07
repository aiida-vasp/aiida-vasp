"""Contains tools to read kpoints and bands from EIGENVALUE files."""

import re

import numpy as np

from .parser import BaseParser


class EigParser(BaseParser):
    """Contains regex and functions to find grammar elements in EIGENVALUE files."""

    def __init__(self, filename):
        res = self.parse_eigenval(filename)
        self.header = res[0]
        self.kpoints = res[1]
        self.bands = res[2]

    @classmethod
    # pylint: disable=too-many-locals
    def parse_eigenval(cls, filename):
        """Parse a VASP EIGENVAL file and extract metadata and a band structure data array"""
        with open(filename) as eig:
            line_0 = cls.line(eig, int)  # read header
            line_1 = cls.line(eig, float)  # "
            line_2 = cls.line(eig, float)  # "
            coord_type = cls.line(eig)  # "
            name = cls.line(eig)  # read name line (can be empty)
            param_0, num_kp, num_bands = cls.line(eig, int)  # read: ? #kp #bands
            data = eig.read()  # rest is data
        num_ions, num_atoms, p00, num_spins = line_0
        data = re.split(cls.empty_line, data)  # list of data blocks
        data = [[line.split() for line in block.splitlines()] for block in data]
        kpoints = np.zeros((num_kp, 4))
        bands = np.zeros((num_spins, num_kp, num_bands))
        for k, field in enumerate(data):  # iterate over data blocks
            kpbs = filter(None, field)  # throw away empty lines
            kpi = map(float, kpbs.pop(0))  # first line of block is kpoints -> pop
            kpoints[k] = kpi
            for point in kpbs:  # rest are band energies
                bands[:, k, int(point[0]) - 1] = point[1:num_spins + 1]  # place energy value in bands[kp, nb] (BandstrucureData format)
        header = {}  # build header dict
        header[0] = line_0
        header[1] = line_1
        header[2] = line_2
        header['n_ions'] = num_ions
        header['n_atoms'] = num_atoms
        header['p00'] = p00
        header['nspin'] = num_spins
        header['cartesian'] = coord_type.startswith(('c', 'C'))
        header['name'] = name
        header['some_num'] = param_0
        header['n_bands'] = num_bands
        header['n_kp'] = num_kp

        return header, kpoints, bands
