"""
EIGENVAL parser.

----------------
The file parser that handles the parsing of EIGENVAL files.
"""

import re

import numpy as np

from aiida_vasp.parsers.file_parsers.parser import BaseFileParser


class EigParser(BaseFileParser):
    """Contains regex and functions to find grammar elements in EIGENVALUE files."""

    PARSABLE_ITEMS = {
        'eigenval-eigenvalues': {
            'inputs': [],
            'name': 'eigenvalues',
            'prerequisites': [],
        },
        'eigenval-kpoints': {
            'inputs': ['structure'],
            'name': 'kpoints',
            'prerequisites': ['structure'],
        },
    }

    def __init__(self, *args, **kwargs):
        super(EigParser, self).__init__(*args, **kwargs)
        self.init_with_kwargs(**kwargs)

    @property
    def _parsed_object(self):
        return self._data_obj

    def _parse_file(self, inputs):
        """Parse a VASP EIGENVAL file and extract metadata and a band structure data array."""

        result = inputs.get('settings', {})
        result = {}

        header, kpoints, bands = self._read_eigenval()
        result['header'] = header
        result['eigenval-eigenvalues'] = bands
        result['eigenval-kpoints'] = kpoints

        return result

    # pylint: disable=too-many-locals
    def _read_eigenval(self):
        """Parse a VASP EIGENVAL file and extract metadata and a band structure data array."""

        with open(self._data_obj.path) as eig:
            line_0 = self.line(eig, int)  # read header
            line_1 = self.line(eig, float)  # "
            line_2 = self.line(eig, float)  # "
            coord_type = self.line(eig)  # "
            name = self.line(eig)  # read name line (can be empty)
            param_0, num_kp, num_bands = self.line(eig, int)  # read: ? #kp #bands
            data = eig.read()  # rest is data
        num_ions, num_atoms, p00, num_spins = line_0
        data = re.split(self.empty_line, data)  # list of data blocks
        data = [[line.split() for line in block.splitlines()] for block in data]
        kpoints = np.zeros((num_kp, 4))
        bands = np.zeros((num_spins, num_kp, num_bands))
        for k, field in enumerate(data):  # iterate over data blocks
            kpbs = [x for x in field if x]  # throw away empty lines
            kpi = [float(x) for x in kpbs.pop(0)]  # first line of block is kpoints -> pop
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
