"""
WAVECAR parser.

---------------
The parser that handles the parsing of WAVECAR.
"""
from aiida_vasp.parsers.object_parsers.parser import BaseFileParser
from aiida_vasp.parsers.node_composer import NodeComposer, get_node_composer_inputs_from_object_parser


class WavecarParser(BaseFileParser):
    """Add WAVECAR as a SinglefileData node."""

    PARSABLE_ITEMS = {
        'wavecar': {
            'inputs': [],
            'name': 'wavecar',
            'prerequisites': []
        },
    }

    def __init__(self, *args, **kwargs):
        super(WavecarParser, self).__init__(*args, **kwargs)
        self._wavecar = None

    def _parse_object(self, inputs):
        """Create a DB Node for the WAVECAR."""
        result = inputs
        result = {}
        wfn = self._data_obj.handler

        if wfn is None:
            return {'wavecar': None}

        result['wavecar'] = wfn
        return result

    @property
    def wavecar(self):
        if self._wavecar is None:
            inputs = get_node_composer_inputs_from_object_parser(self)
            self._wavecar = NodeComposer.compose('vasp.wavefun', inputs)
        return self._wavecar
