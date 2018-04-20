"""Tools for parsing CHGCAR files."""

from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.io.parser import BaseFileParser


class WavecarParser(BaseFileParser):
    """Add WAVECAR as a single file node."""

    PARSABLE_ITEMS = {
        'wavecar': {
            'inputs': [],
            'parsers': ['WAVECAR'],
            'nodeName': 'wavecar',
            'prerequisites': []
        },
    }

    def __init__(self, path, filename, cls):
        super(WavecarParser, self).__init__(cls)
        self._filepath = path
        self._filename = filename
        self._parsable_items = WavecarParser.PARSABLE_ITEMS
        self._parsed_data = {}

    def _parse_file(self, inputs):
        """Create a DB Node for the WAVECAR file"""
        result = inputs
        result = {}
        wfn = self._filepath

        if wfn is None:
            return {'wavecar': None}

        result['wavecar'] = get_data_class('vasp.wavefun')()
        result['wavecar'].set_file(wfn)
        return result
