"""Tools for parsing CHGCAR files."""

from aiida.orm import DataFactory
from aiida_vasp.io.parser import BaseParser


class WavecarParser(BaseParser):
    """Add WAVECAR as a single file node."""

    PARSABLE_ITEMS = {
        'wavecar': {
            'inputs': [],
            'parsers': ['WAVECAR'],
            'nodeName': 'wavecar',
            'prerequisites': []
        },
    }

    def __init__(self, path, filename):
        super(WavecarParser, self).__init__()
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

        result['wavecar'] = DataFactory('vasp.wavefun')()
        result['wavecar'].set_file(wfn)
        return result
