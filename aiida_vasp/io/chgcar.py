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
            'chgcar': None,
        },
    }

    def __init__(self, path, filename):
        super(ChgcarParser, self).__init__()
        self._filepath = path
        self._filename = filename
        self._parsable_items = ChgcarParser.PARSABLE_ITEMS

    def _get_chgcar(self, inputs):
        """Create a DB Node for the CHGCAR file"""
        result = inputs

        chgcar = self._filepath
        if chgcar is None:
            return {'chgcar': None}

        result = DataFactory('vasp.chargedensity')()
        result.set_file(chgcar)

        return {'chgcar': result}
