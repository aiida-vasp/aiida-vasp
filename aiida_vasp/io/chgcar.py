"""Tools for parsing CHGCAR files."""

from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.io.parser import BaseFileParser


class ChgcarParser(BaseFileParser):
    """Add CHGCAR as a single file node."""

    PARSABLE_ITEMS = {
        'chgcar': {
            'inputs': [],
            'nodeName': 'chgcar',
            'prerequisites': []
        },
    }

    def __init__(self, *args, **kwargs):
        super(ChgcarParser, self).__init__(*args, **kwargs)
        self.init_with_kwargs(**kwargs)

    def _parse_file(self, inputs):
        """Create a DB Node for the CHGCAR file"""
        result = inputs
        result = {}

        chgcar = self._data_obj.path
        if chgcar is None:
            return {'chgcar': None}

        result['chgcar'] = get_data_class('vasp.chargedensity')()
        result['chgcar'].set_file(chgcar)

        return result
