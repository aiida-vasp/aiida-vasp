"""
WAVECAR parser.

---------------
The file parser that handles the parsing of WAVECAR files.
"""
from aiida_vasp.parsers.file_parsers.parser import BaseFileParser
from aiida_vasp.parsers.node_composer import NodeComposer, get_node_composer_inputs_from_file_parser


class WavecarParser(BaseFileParser):
    """Add WAVECAR as a single file node."""

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

    def _parse_file(self, inputs):
        """Create a DB Node for the WAVECAR file."""
        result = inputs
        result = {}
        wfn = self._data_obj.path

        if wfn is None:
            return {'wavecar': None}

        result['wavecar'] = wfn
        return result

    @property
    def wavecar(self):
        if self._wavecar is None:
            inputs = get_node_composer_inputs_from_file_parser(self)
            self._wavecar = NodeComposer.compose('vasp.wavefun', inputs)
        return self._wavecar
