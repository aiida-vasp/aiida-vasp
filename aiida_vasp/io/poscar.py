"""
Tools for parsing POSCAR files.
"""
from aiida.orm import DataFactory
from aiida_vasp.io.parser import BaseParser


class PoscarParser(BaseParser):
    """
    Parse a POSCAR format file into a StructureData node.
    """

    PARSABLE_ITEMS = {
        'structure': {
            'structure': None,
        },
    }

    def __init__(self, path, filename):
        super(PoscarParser, self).__init__()
        self._filepath = path
        self._filename = filename
        self._parsable_items = PoscarParser.PARSABLE_ITEMS

    def _get_structure(self, inputs):
        '''read POSCAR format file for output structure'''
        from ase.io import read

        result = inputs
        result = DataFactory('structure')()
        cont = self._filepath
        if not cont:
            return {'structure': None}
        result.set_ase(read(cont, format='vasp'))

        return {'structure': result}
