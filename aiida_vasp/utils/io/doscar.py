"""DOSCAR (VASP format) utilities"""
import numpy as np

from .parser import BaseParser


class DosParser(BaseParser):
    """
    parse a DOSCAR file from a vasp run
    """

    def __init__(self, filename, **kwargs):
        self.ispin = kwargs.get('ispin')
        self.lorbit = kwargs.get('lorbit')
        self.rwigs = kwargs.get('rwigs')
        self.header, self.tdos, self.pdos = self.parse_doscar(filename)

    @classmethod
    # pylint: disable=too-many-locals
    def parse_doscar(cls, filename):
        """Read a VASP DOSCAR file and extract metadata and a density of states data array"""
        with open(filename) as dos:
            num_ions, num_atoms, p00, p01 = cls.line(dos, int)
            line_0 = cls.line(dos, float)
            line_1 = cls.line(dos, float)
            coord_type = cls.line(dos)
            sys = cls.line(dos)
            line_2 = cls.line(dos, float)
            emax, emin, ndos, efermi, weight = line_2
            ndos = int(ndos)
            raw = cls.splitlines(dos)

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

        return header, tdos, pdos
