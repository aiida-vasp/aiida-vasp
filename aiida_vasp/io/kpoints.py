"""Utils for VASP KPOINTS format"""
import numpy as np

from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.io.parser import BaseFileParser


class KpParser(BaseFileParser):
    """
    Parser for VASP KPOINTS format

    This is a remainder of the previous VaspParser and at the moment only capable of
    parsing KPOINTS files containing an explicit list of k-points. This should be
    replaced by a more capable parser in the future.
    """

    PARSABLE_ITEMS = {
        'kpoints': {
            'inputs': [],
            'parsers': ['EIGENVAL', 'IBZKPT'],
            'nodeName': 'kpoints',
            'prerequisites': []
        },
    }

    def __init__(self, *args, **kwargs):
        super(KpParser, self).__init__(*args, **kwargs)
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
        header['cartesian'] = fobj_or_str.readline().startswith(('c', 'C', 'k', 'K'))
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

        header, res = self.parse_kp_file(self._file_path)
        kpoints = res[:, :3]
        if res.shape[1] == 4:
            weights = res[:, 3]
        else:
            weights = None

        kpout = get_data_class('array.kpoints')()
        kpout.set_kpoints(kpoints, weights=weights, cartesian=header['cartesian'])

        result['kpoints_header'] = header
        result['kpoints_raw'] = kpoints
        result['weights'] = weights
        result['kpoints'] = kpout

        return result
