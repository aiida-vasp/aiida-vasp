"""
DOSCAR parser.

--------------
The file parser that handles the parsing of DOSCAR files.
"""
# pylint: disable=unsubscriptable-object  # pylint/issues/3139
import numpy as np

from aiida_vasp.parsers.node_composer import NodeComposer, get_node_composer_inputs_from_file_parser
from aiida_vasp.parsers.file_parsers.parser import BaseFileParser

# Map from number of columns in DOSCAR to dtype.
DTYPES_TDOS = {
    3: np.dtype([('energy', float), ('total', float), ('integrated', float)]),
    5: np.dtype([('energy', float), ('total', float, (2,)), ('integrated', float, (2,))]),
}

DTYPES_PDOS = {
    # l-decomposed
    4:
        np.dtype([('energy', float), ('s', float), ('p', float), ('d', float)]),
    7:
        np.dtype([('energy', float), ('s', float, (2,)), ('p', float, (2,)), ('d', float, (2,))]),
    13:
        np.dtype([('energy', float), ('s', float, (4,)), ('p', float, (4,)), ('d', float, (4,))]),
    5:
        np.dtype([('energy', float), ('s', float), ('p', float), ('d', float), ('f', float)]),
    9:
        np.dtype([('energy', float), ('s', float, (2,)), ('p', float, (2,)), ('d', float, (2,)), ('f', float, (2,))]),
    17:
        np.dtype([('energy', float), ('s', float, (4,)), ('p', float, (4,)), ('d', float, (4,)), ('f', float, (4,))]),
    # lm-decomposed
    10:
        np.dtype([('energy', float), ('s', float), ('py', float), ('px', float), ('pz', float), ('dxy', float), ('dyz', float),
                  ('dz2', float), ('dxz', float), ('dx2-y2', float)]),
    18:
        np.dtype([('energy', float), ('s', float), ('py', float), ('px', float), ('pz', float), ('dxy', float), ('dyz', float),
                  ('dz2', float), ('dxz', float), ('dx2-y2', float), ('fy(3x2-y2)', float), ('fxyz', float), ('fyz2', float),
                  ('fz3', float), ('fxz2', float), ('fz(x2-y2)', float), ('fx(x2-3y2)', float)]),
    19:
        np.dtype([('energy', float), ('s', float, (2,)), ('py', float, (2,)), ('px', float, (2,)), ('pz', float, (2,)),
                  ('dxy', float, (2,)), ('dyz', float, (2,)), ('dz2', float, (2,)), ('dxz', float, (2,)), ('dx2-y2', float, (2,))]),
    35:
        np.dtype([('energy', float), ('s', float, (2,)), ('py', float, (2,)), ('px', float, (2,)), ('pz', float, (2,)),
                  ('dxy', float, (2,)), ('dyz', float, (2,)), ('dz2', float, (2,)), ('dxz', float, (2,)), ('dx2-y2', float, (2,)),
                  ('fy(3x2-y2)', float, (2,)), ('fxyz', float, (2,)), ('fyz2', float, (2,)), ('fz3', float, (2,)), ('fxz2', float, (2,)),
                  ('fz(x2-y2)', float, (2,)), ('fx(x2-3y2)', float, (2,))]),
    37:
        np.dtype([('energy', float), ('s', float, (4,)), ('py', float, (4,)), ('px', float, (4,)), ('pz', float, (4,)),
                  ('dxy', float, (4,)), ('dyz', float, (4,)), ('dz2', float, (4,)), ('dxz', float, (4,)), ('x2-y2', float, (4,))]),
    69:
        np.dtype([('energy', float), ('s', float, (4,)), ('py', float, (4,)), ('px', float, (4,)), ('pz', float, (4,)),
                  ('dxy', float, (4,)), ('dyz', float, (4,)), ('dz2', float, (4,)), ('dxz', float, (4,)), ('dx2-y2', float, (4,)),
                  ('fy(3x2-y2)', float, (4,)), ('fxyz', float, (4,)), ('fyz2', float, (4,)), ('fz3', float, (4,)), ('fxz2', float, (4,)),
                  ('fz(x2-y2)', float, (4,)), ('fx(x2-3y2)', float, (4,))]),
}

# Mapping between the number of columns to the number of spins
COLSPIN_MAP = {7: 2, 9: 2, 19: 2, 35: 2, 13: 4, 17: 4, 37: 4, 69: 4, 4: 1, 5: 1, 10: 1, 18: 1}


class DosParser(BaseFileParser):
    """Parse a DOSCAR file from a vasp run."""

    PARSABLE_ITEMS = {
        'doscar-dos': {
            'inputs': [],
            'name': 'dos',
            'prerequisites': [],
        },
    }

    def __init__(self, *args, **kwargs):
        super(DosParser, self).__init__(*args, **kwargs)
        self._dos = None

    def _parse_file(self, inputs):
        """Read a VASP DOSCAR file and extract metadata and a density of states data array."""

        result = inputs
        result = {}

        header, pdos, tdos = self._read_doscar()

        result['doscar-dos'] = {}
        result['header'] = header
        if pdos.size != 0:
            result['doscar-dos']['pdos'] = pdos
        if tdos.size != 0:
            result['doscar-dos']['tdos'] = tdos

        # No data is avaliable - return None
        if not result['doscar-dos']:
            return {'doscar-dos': None}

        return result

    # pylint: disable=too-many-locals, too-many-statements
    def _read_doscar(self):
        """Read a VASP DOSCAR file and extract metadata and a density of states data array."""

        with open(self._data_obj.path) as dos:
            num_ions, num_atoms, p00, p01 = self.line(dos, int)
            line_0 = self.line(dos, float)
            line_1 = self.line(dos, float)
            coord_type = self.line(dos)
            sys = self.line(dos)
            line_2 = self.line(dos, float)
            emax, emin, ndos, efermi, weight = line_2
            ndos = int(ndos)
            raw = self.splitlines(dos)

        # Get the number of columns for the tdos section.
        count = len(raw[ndos - 1])

        num_spin = 1
        if count == 5:
            num_spin = 2

        tdos_raw = np.array(raw[:ndos])
        tdos = np.zeros((tdos_raw.shape[0]), DTYPES_TDOS[count])
        tdos['energy'] = tdos_raw[:, 0]
        for i, name in enumerate(DTYPES_TDOS[count].names[1:]):
            if num_spin == 1:
                tdos[name] = np.squeeze(tdos_raw[:, i + 1:i + 1 + num_spin], axis=1)
            else:
                tdos[name] = tdos_raw[:, i + 1:i + 1 + num_spin]

        pdos_items = []
        pdos = np.array([])  # Empty array by default
        if line_2 in raw:
            for _ in range(num_ions):
                start = raw.index(line_2) + 1
                pdos_items += [raw[start:start + ndos]]

            # Get the number of columns for the pdos section.
            count = len(pdos_items[-1][-1])
            pdos_raw = np.array(pdos_items)

            # Adjust the spin according to the column definitions
            num_spin = COLSPIN_MAP.get(count)
            if num_spin is None:
                raise ValueError(f'Unkown column count: {count} in DOSCAR')

            pdos = np.zeros((pdos_raw.shape[0], pdos_raw.shape[1]), DTYPES_PDOS[count])
            pdos['energy'] = pdos_raw[:, :, 0]
            for i, name in enumerate(DTYPES_PDOS[count].names[1:]):
                if num_spin == 1:  # Only squeeze if there is only one spin component
                    pdos[name] = np.squeeze(pdos_raw[:, :, i + 1:i + 1 + num_spin], axis=2)
                else:
                    pdos[name] = pdos_raw[:, :, i + 1:i + 1 + num_spin]

        header = {}
        header[0] = line_0
        header[1] = line_1
        header[2] = line_2
        header['n_ions'] = num_ions
        header['n_atoms'] = num_atoms
        header['p00'] = p00
        header['p01'] = p01
        header['cartesian'] = coord_type.startswith(('c', 'C'))
        header['name'] = sys
        header['emax'] = emax
        header['emin'] = emin
        header['n_dos'] = ndos
        header['efermi'] = efermi
        header['weight'] = weight

        return header, pdos, tdos

    @property
    def dos(self):
        if self._dos is None:
            inputs = get_node_composer_inputs_from_file_parser(self, quantity_keys=['doscar-dos'])
            self._dos = NodeComposer.compose('array', inputs)
        return self._dos
