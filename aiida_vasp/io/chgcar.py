"""Tools for parsing CHGCAR files."""

from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.io.parser import BaseFileParser


class ChgcarParser(BaseFileParser):
    """Add CHGCAR as a single file node."""

    PARSABLE_ITEMS = {
        'chgcar': {
            'inputs': [],
            'parsers': ['CHGCAR'],
            'nodeName': 'chgcar',
            'prerequisites': []
        },
    }

    def __init__(self, path, filename, cls):
        super(ChgcarParser, self).__init__(cls)
        self._filepath = path
        self._filename = filename
        self._parsable_items = ChgcarParser.PARSABLE_ITEMS
        self._parsed_data = {}

    def _parse_file(self, inputs):
        """Create a DB Node for the CHGCAR file"""
        result = inputs
        result = {}

        chgcar = self._filepath
        if chgcar is None:
            return {'chgcar': None}

        result['chgcar'] = get_data_class('vasp.chargedensity')()
        result['chgcar'].set_file(chgcar)

        return result
