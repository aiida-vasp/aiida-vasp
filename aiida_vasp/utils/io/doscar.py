from parser import BaseParser
import numpy as np


class DosParser(BaseParser):
    '''
    parse a DOSCAR file from a vasp run
    '''
    def __init__(self, filename, **kwargs):
        self.ispin = kwargs.get('ispin')
        self.lorbit = kwargs.get('lorbit')
        self.rwigs = kwargs.get('rwigs')
        h, t, p = self.parse_doscar(filename)
        self.header = h
        self.tdos = t
        self.pdos = p

    @classmethod
    def parse_doscar(cls, filename):
        with open(filename) as dos:
            ni, na, p00, p01 = cls.line(dos, int)
            l0 = cls.line(dos, float)
            l1 = cls.line(dos, float)
            coord_type = cls.line(dos)
            sys = cls.line(dos)
            l2 = cls.line(dos, float)
            emax, emin, ndos, efermi, weight = l2
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
        if l2 in raw:
            for i in range(ni):
                start = raw.index(l2) + 1
                pdos += [raw[start:start+ndos]]
            pdos = np.array(pdos)

        header = {}
        header[0] = l0
        header[1] = l1
        header[2] = l2
        header['n_ions'] = ni
        header['n_atoms'] = na
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
