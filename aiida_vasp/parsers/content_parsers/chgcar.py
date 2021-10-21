"""
CHGCAR parser.

--------------
The object parser that handles the parsing of CHGCAR.
"""

from aiida_vasp.parsers.content_parsers.parser import BaseFileParser
from aiida_vasp.parsers.node_composer import NodeComposer, get_node_composer_inputs_from_object_parser


class ChgcarParser(BaseFileParser):
    """Add CHGCAR as a SinglefileData node."""

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

    def _parse_object(self, inputs):
        """Create a DB Node for CHGCAR."""
        result = inputs
        result = {}

        chgcar = self._data_obj.handler
        if chgcar is None:
            return {'chgcar': None}

        result['chgcar'] = chgcar

        return result

    @property
    def chgcar(self):
        if self._chgcar is None:
            inputs = get_node_composer_inputs_from_object_parser(object_parser=self)
            self._chgcar = NodeComposer.compose('vasp.chargedensity', inputs)
        return self._chgcar
