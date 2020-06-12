"""
WAVECAR parser.

---------------
The file parser that handles the parsing of WAVECAR files.
"""
from aiida_vasp.parsers.file_parsers.parser import BaseFileParser
from aiida_vasp.parsers.node_composer import NodeComposer


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
        self.init_with_kwargs(**kwargs)

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
            composer = NodeComposer(file_parsers=[self])
            self._wavecar = composer.compose('vasp.wavefun')
        return self._wavecar
