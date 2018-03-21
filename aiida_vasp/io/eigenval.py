"""Contains tools to read kpoints and bands from EIGENVALUE files."""

import re

import numpy as np

from aiida.orm import DataFactory
from aiida_vasp.io.parser import BaseParser


class EigParser(BaseParser):
    """Contains regex and functions to find grammar elements in EIGENVALUE files."""

    PARSABLE_ITEMS = {
        'bands': {
            'inputs': ['structure', 'kpoints', 'occupations'],
            'parsers': ['EIGENVAL', 'vasprun.xml'],
            'nodeName': 'bands',
            'prerequisites': ['structure', 'occupations']
        },
    }

    def __init__(self, path, filename):
        super(EigParser, self).__init__()
        self._filepath = path
        self._filename = filename
        self._parsed_data = None
        self._parsable_items = EigParser.PARSABLE_ITEMS

    def _parse_file(self, inputs):
        """Parse a VASP EIGENVAL file and extract metadata and a band structure data array"""

        result = inputs.get('settings', {})
        result = {}

        header, kpoints, bands = self._read_eigenval()

        result['header'] = header

        bsnode = DataFactory('array.bands')()
        kpout = DataFactory('array.kpoints')()

        structure = inputs.get('structure')
        if structure is None:
            return {'bands': None, 'kpoints': None}

        bsnode.set_cell(structure.get_ase().get_cell())
        kpout.set_cell(structure.get_ase().get_cell())

        kpoints_inp = inputs.get('kpoints')
        if kpoints_inp:
            bsnode.set_kpointsdata(kpoints_inp)

            if kpoints_inp.labels:
                bsnode.labels = kpoints_inp.labels
        else:
            bsnode.set_kpoints(kpoints[:, :3], weights=kpoints[:, 3], cartesian=False)

        bsnode.set_bands(bands, occupations=inputs['occupations'])
        kpout.set_kpoints(kpoints[:, :3], weights=kpoints[:, 3], cartesian=False)

        result['bands'] = bsnode
        result['kpoints'] = kpout

        return result

    # pylint: disable=too-many-locals
    def _read_eigenval(self):
        """Parse a VASP EIGENVAL file and extract metadata and a band structure data array"""

        with open(self._filepath) as eig:
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
