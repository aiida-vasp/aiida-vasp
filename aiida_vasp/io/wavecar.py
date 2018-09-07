"""Tools for parsing CHGCAR files."""

from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.io.parser import BaseFileParser


class WavecarParser(BaseFileParser):
    """Add WAVECAR as a single file node."""

    PARSABLE_ITEMS = {
        'wavecar': {
            'inputs': [],
            'nodeName': 'wavecar',
            'prerequisites': []
        },
    }

    def __init__(self, *args, **kwargs):
        super(WavecarParser, self).__init__(*args, **kwargs)
        self.init_with_kwargs(**kwargs)

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
