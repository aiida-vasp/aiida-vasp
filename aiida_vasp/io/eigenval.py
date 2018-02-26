"""
This module contains tools to read kpoints and bands from EIGENVALUE files.
"""

import re

import numpy as np

from aiida.orm import DataFactory
from aiida_vasp.io.parser import BaseParser


class EigParser(BaseParser):
    """
    contains regex and functions to find grammar elements
    in EIGENVALUE files
    """

    PARSABLE_ITEMS = {
        'bands': {
            'inputs': ['structure', 'kpoints', 'occupations'],
            'parsers': ['EIGENVAL', 'vasprun.xml'],
            'nodeName': 'bands',
            'prerequesites': ['structure', 'occupations']
        },
    }

    def __init__(self, path, filename):
        super(EigParser, self).__init__()
        self._filepath = path
        self._filename = filename
        self._parsed_data = None
        self._parsable_items = EigParser.PARSABLE_ITEMS

    def _get_bands(self, inputs):
        '''
        Create a bands and a kpoints node from values in eigenvalue.

        returns: bsnode, kpout
        - bsnode: BandsData containing eigenvalues from EIGENVAL
                and occupations from vasprun.xml
        - kpout: KpointsData containing kpoints from EIGENVAL,

        both bsnode as well as kpnode come with cell unset
        '''
        eig = self._filepath
        if not eig:
            return {'bands': None, 'kpoints': None}

        _, kpoints, bands = self.parse_eigenval(eig)

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

        return {'bands': bsnode, 'kpoints': kpout}

    def parse_eigenval(self, filename):
        """Return values from EIGENVAL. Either return cached values or parse the file."""

        if not self._parsed_data:
            self._parsed_data = {}
            self._parsed_data['header'], self._parsed_data['kpoints'], self._parsed_data['bands'] = self._parse_eigenval(filename)

        return self._parsed_data['header'], self._parsed_data['kpoints'], self._parsed_data['bands']

    # pylint: disable=too-many-locals
    def _parse_eigenval(self, filename):
        """Parse a VASP EIGENVAL file and extract metadata and a band structure data array"""

        with open(filename) as eig:
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
