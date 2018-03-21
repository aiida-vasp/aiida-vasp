"""Utils for VASP KPOINTS format"""
import numpy as np

from aiida.orm import DataFactory
from aiida_vasp.io.parser import BaseParser


class KpParser(BaseParser):
    """Parser for VASP KPOINTS format"""

    PARSABLE_ITEMS = {
        'kpoints': {
            'inputs': [],
            'parsers': ['EIGENVAL', 'IBZKPT'],
            'nodeName': 'kpoints',
            'prerequisites': []
        },
    }

    def __init__(self, path, filename):
        super(KpParser, self).__init__()
        self._filepath = path
        self._filename = filename
        self._parsable_items = KpParser.PARSABLE_ITEMS
        self._parsed_data = {}

    @classmethod
    def parse_kp(cls, fobj_or_str):
        """Parse VASP KPOINTS files"""

        if isinstance(fobj_or_str, str):
            from StringIO import StringIO
            fobj_or_str = StringIO(fobj_or_str)
        header = {}
        header['name'] = fobj_or_str.readline()
        header['nkp'] = cls.line(fobj_or_str, d_type=int)
        header['cartesian'] = not fobj_or_str.readline().startswith(('r', 'R'))
        lines = np.array(cls.splitlines(fobj_or_str))
        return header, lines

    @classmethod
    def parse_kp_file(cls, fname):
        with open(fname) as kpoints_f:
            return cls.parse_kp(kpoints_f)

    def _parse_file(self, inputs):
        """Create a DB Node for the IBZKPT file"""

        result = inputs
        result = {}

        header, res = self.parse_kp_file(self._filepath)
        kpoints = res[:, :3]
        if res.shape[1] == 4:
            weights = res[:, 3]
        else:
            weights = None

        kpout = DataFactory('array.kpoints')()
        kpout.set_kpoints(kpoints, weights=weights, cartesian=header['cartesian'])

        result['kpoints'] = kpout

        return result
