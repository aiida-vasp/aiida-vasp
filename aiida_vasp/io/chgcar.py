"""
Tools for parsing CHGCAR files.
"""
from aiida.orm import DataFactory
from aiida_vasp.io.parser import BaseParser


class ChgcarParser(BaseParser):
    """
    Add CHGCAR as a single file node.
    """

    PARSABLE_ITEMS = {
        'chgcar': {
            'inputs': [],
            'parsers': ['CHGCAR'],
            'nodeName': 'chgcar',
            'prerequisites': []
        },
    }

    def __init__(self, path, filename):
        super(ChgcarParser, self).__init__()
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

        result['chgcar'] = DataFactory('vasp.chargedensity')()
        result['chgcar'].set_file(chgcar)

        return result
