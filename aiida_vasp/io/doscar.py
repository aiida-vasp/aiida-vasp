"""DOSCAR (VASP format) utilities"""
import numpy as np

from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.io.parser import BaseFileParser

# Map from number of columns in DOSCAR to dtype.
DTYPES = {
    3:
    np.dtype([('energy', float), ('total', float), ('integrated', float)]),
    5:
    np.dtype([('energy', float), ('total', float, (2,)), ('integrated', float, (2,))]),
    10:
    np.dtype([('energy', float), ('s', float), ('py', float), ('px', float), ('pz', float), ('dxy', float), ('dyz', float), ('dz2', float),
              ('dxz', float), ('x2-y2', float)]),
    19:
    np.dtype([('energy', float), ('s', float, (2,)), ('py', float, (2,)), ('px', float, (2,)), ('pz', float, (2,)), ('dxy', float, (2,)),
              ('dyz', float, (2,)), ('dz2', float, (2,)), ('dxz', float, (2,)), ('x2-y2', float, (2,))]),
    37:
    np.dtype([('energy', float), ('s', float, (4,)), ('py', float, (4,)), ('px', float, (4,)), ('pz', float, (4,)), ('dxy', float, (4,)),
              ('dyz', float, (4,)), ('dz2', float, (4,)), ('dxz', float, (4,)), ('x2-y2', float, (4,))])
}


class DosParser(BaseFileParser):
    """Parse a DOSCAR file from a vasp run."""

    PARSABLE_ITEMS = {
        'doscar-dos': {
            'inputs': [],
            'nodeName': 'dos',
            'prerequisites': [],
            'alternatives': ['dos']
        },
    }

    def __init__(self, *args, **kwargs):
        super(DosParser, self).__init__(*args, **kwargs)
        self.init_with_kwargs(**kwargs)

    def _parse_file(self, inputs):
        """Read a VASP DOSCAR file and extract metadata and a density of states data array"""

        result = inputs
        result = {}

        header, pdos, tdos = self._read_doscar()

        result['header'] = header

        for array in [pdos, tdos]:
            if array.size == 0:
                return {'doscar-dos': None}

        dosnode = get_data_class('array')()
        dosnode.set_array('pdos', pdos)
        dosnode.set_array('tdos', tdos)

        result['doscar-dos'] = dosnode

        return result

    # pylint: disable=too-many-locals
    def _read_doscar(self):
        """Read a VASP DOSCAR file and extract metadata and a density of states data array"""

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
        if count == 9:
            num_spin = 4

        tdos_raw = np.array(raw[:ndos])
        tdos = np.zeros((tdos_raw.shape[0]), DTYPES[count])
        tdos['energy'] = tdos_raw[:, 0]
        for i, name in enumerate(DTYPES[count].names[1:]):
            tdos[name] = np.squeeze(tdos_raw[:, i + 1:i + 1 + num_spin], axis=1)

        pdos = []
        if line_2 in raw:
            for _ in range(num_ions):
                start = raw.index(line_2) + 1
                pdos += [raw[start:start + ndos]]

            # Get the number of columns for the pdos section.
            count = len(pdos[-1][-1])
            pdos_raw = np.array(pdos)

            pdos = np.zeros((pdos_raw.shape[0], pdos_raw.shape[1]), DTYPES[count])
            pdos['energy'] = pdos_raw[:, :, 0]
            for i, name in enumerate(DTYPES[count].names[1:]):
                pdos[name] = np.squeeze(pdos_raw[:, :, i + 1:i + 1 + num_spin], axis=2)

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
