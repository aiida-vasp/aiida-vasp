"""
CHGCAR parser.

--------------
The file parser that handles the parsing of CHGCAR files.
"""

from aiida_vasp.parsers.file_parsers.parser import BaseFileParser
from aiida_vasp.parsers.node_composer import NodeComposer


class ChgcarParser(BaseFileParser):
    """Add CHGCAR as a single file node."""

    PARSABLE_ITEMS = {
        'chgcar': {
            'inputs': [],
            'name': 'chgcar',
            'prerequisites': []
        },
    }

    def __init__(self, *args, **kwargs):
        super(ChgcarParser, self).__init__(*args, **kwargs)
        self._chgcar = None
        self.init_with_kwargs(**kwargs)

    def _parse_file(self, inputs):
        """Create a DB Node for the CHGCAR file."""
        result = inputs
        result = {}

        chgcar = self._data_obj.path
        if chgcar is None:
            return {'chgcar': None}

        result['chgcar'] = chgcar

        return result

    @property
    def chgcar(self):
        if self._chgcar is None:
            composer = NodeComposer(file_parsers=[self])
            self._chgcar = composer.compose('vasp.chargedensity')
        return self._chgcar
