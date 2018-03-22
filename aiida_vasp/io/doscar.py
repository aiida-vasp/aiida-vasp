"""DOSCAR (VASP format) utilities"""
import numpy as np

from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.io.parser import BaseParser


class DosParser(BaseParser):
    """Parse a DOSCAR file from a vasp run."""

    PARSABLE_ITEMS = {
        'dos': {
            'inputs': ['vrp_pdos', 'vrp_tdos'],
            'parsers': ['vasprun.xml', 'DOSCAR'],
            'nodeName': 'dos',
            'prerequisites': ['vrp_pdos', 'vrp_tdos']
        },
    }

    def __init__(self, path, filename):
        super(DosParser, self).__init__()
        self._filepath = path
        self._filename = filename
        self._parsable_items = DosParser.PARSABLE_ITEMS
        self._parsed_data = {}

    def _parse_file(self, inputs):
        """Read a VASP DOSCAR file and extract metadata and a density of states data array"""

        result = inputs
        result = {}

        header, dcp_pdos, dcp_tdos = self._read_doscar()

        result['header'] = header

        vrp_pdos = inputs.get('vrp_pdos', np.array([]))
        vrp_tdos = inputs.get('vrp_tdos', np.array([]))

        for array in [vrp_pdos, vrp_tdos, dcp_pdos, dcp_tdos]:
            if array.size == 0:
                return {'dos': None}

        dosnode = get_data_class('array')()
        # vrp_pdos is a numpy array, and thus not directly bool-convertible
        if vrp_pdos.size > 0:
            pdos = vrp_pdos.copy()
            for i, name in enumerate(vrp_pdos.dtype.names[1:]):
                num_spins = vrp_pdos.shape[1]
                # ~ pdos[name] = dcp[:, :, i+1:i+1+ns].transpose(0,2,1)
                cur = dcp_pdos[:, :, i + 1:i + 1 + num_spins].transpose(0, 2, 1)
                cond = vrp_pdos[name] < 0.1
                pdos[name] = np.where(cond, cur, vrp_pdos[name])
            dosnode.set_array('pdos', pdos)
        num_spins = 1
        if dcp_tdos.shape[1] == 5:
            num_spins = 2
        tdos = vrp_tdos[:num_spins, :].copy()
        for i, name in enumerate(vrp_tdos.dtype.names[1:]):
            cur = dcp_tdos[:, i + 1:i + 1 + num_spins].transpose()
            cond = vrp_tdos[:num_spins, :][name] < 0.1
            tdos[name] = np.where(cond, cur, vrp_tdos[:num_spins, :][name])
        dosnode.set_array('tdos', tdos)
        result['dos'] = dosnode

        return result

    # pylint: disable=too-many-locals
    def _read_doscar(self):
        """Read a VASP DOSCAR file and extract metadata and a density of states data array"""

        with open(self._filepath) as dos:
            num_ions, num_atoms, p00, p01 = self.line(dos, int)
            line_0 = self.line(dos, float)
            line_1 = self.line(dos, float)
            coord_type = self.line(dos)
            sys = self.line(dos)
            line_2 = self.line(dos, float)
            emax, emin, ndos, efermi, weight = line_2
            ndos = int(ndos)
            raw = self.splitlines(dos)

        # either (e tot intd) or (e tot^ tot_ intd^ intd_)
        tdos = np.array(raw[:ndos])
        # either (e s (p) (d)) -> 10
        # or (e s^ s_ (p^ p_) (d^ d_)) -> 19
        # or (e s[m] (p[m]) (d[m]) -> 37
        # or (e s^[m] s_[m] (p^[m] p_[m]) (d^[m] d_[m]) -> 73
        # probably format later with vasprun or PROCAR info?
        # from vasprun: pdos[i][1+j::n_spin] <-> vrunpdos[i][j][1:]
        pdos = []
        if line_2 in raw:
            for _ in range(num_ions):
                start = raw.index(line_2) + 1
                pdos += [raw[start:start + ndos]]
            pdos = np.array(pdos)

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
