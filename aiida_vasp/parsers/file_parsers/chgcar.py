"""
CHGCAR parser.

--------------
The file parser that handles the parsing of CHGCAR files.
"""

from aiida_vasp.parsers.file_parsers.parser import BaseFileParser
from aiida_vasp.parsers.node_composer import NodeComposer, get_node_composer_inputs


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
            inputs = get_node_composer_inputs(node_type='vasp.chargedensity', file_parser=self)
            self._chgcar = NodeComposer.compose('vasp.chargedensity', inputs)
        return self._chgcar
