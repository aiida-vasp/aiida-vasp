"""
Tools for parsing CHGCAR files.
"""
from aiida.orm import DataFactory
from aiida_vasp.io.parser import BaseParser


class WavecarParser(BaseParser):
    """
    Add CHGCAR as a single file node.
    """

    PARSABLE_ITEMS = {
        'wavecar': {
            'inputs': [],
            'parsers': ['WAVECAR'],
            'nodeName': 'wavecar',
            'prerequesites': []
        },
    }

    def __init__(self, path, filename):
        super(WavecarParser, self).__init__()
        self._filepath = path
        self._filename = filename
        self._parsable_items = WavecarParser.PARSABLE_ITEMS

    def _get_wavecar(self, inputs):
        """Create a DB Node for the WAVECAR file"""
        result = inputs
        wfn = self._filepath

        if wfn is None:
            return {'wavecar': None}

        result = DataFactory('vasp.wavefun')()
        result.set_file(wfn)
        return {'wavecar': result}
