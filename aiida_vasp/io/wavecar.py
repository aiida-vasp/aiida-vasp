"""Tools for parsing CHGCAR files."""

from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.io.parser import BaseFileParser, SingleFile


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

    def __init__(self, *args, **kwargs):
        super(WavecarParser, self).__init__(*args, **kwargs)
        self.init_with_kwargs(**kwargs)

    def _init_with_filepath(self, filepath):
        self._data_obj = SingleFile(filepath=filepath)
        self._parsable_items = WavecarParser.PARSABLE_ITEMS
        self._parsed_data = {}

    def _init_with_data(self, data):
        self._data_obj = SingleFile(data=data)
        self._parsed_data = {}

    @property
    def _parsed_object(self):
        return self._data_obj

    def _parse_file(self, inputs):
        """Create a DB Node for the WAVECAR file"""
        result = inputs
        result = {}
        wfn = self._data_obj.path

        if wfn is None:
            return {'wavecar': None}

        result['wavecar'] = get_data_class('vasp.wavefun')()
        result['wavecar'].set_file(wfn)
        return result
