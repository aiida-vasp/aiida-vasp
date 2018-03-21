"""Tools for parsing OUTCAR files."""

from aiida.orm import DataFactory
from aiida_vasp.io.parser import BaseParser

DEFAULT_OPTIONS = {'quantities_to_parse': ['volume', 'energies', 'efermi']}


class OutcarParser(BaseParser):
    """Parse OUTCAR into a dictionary, which is supposed to be turned into ParameterData later."""

    FILE_NAME = 'OUTCAR'
    PARSABLE_ITEMS = {
        'volume': {
            'inputs': ['parameters'],
            'parsers': ['OUTCAR'],
            'nodeName': 'parameters',
            'prerequisites': []
        },
        'energies': {
            'inputs': ['parameters'],
            'parsers': ['OUTCAR'],
            'nodeName': 'parameters',
            'prerequisites': []
        },
        'efermi': {
            'inputs': ['parameters'],
            'parsers': ['OUTCAR'],
            'nodeName': 'parameters',
            'prerequisites': []
        },
        'parameters': {
            'inputs': [],
            'parsers': ['OUTCAR'],
            'nodeName': 'parameters',
            'prerequisites': []
        },
    }

    def __init__(self, path, filename):
        super(OutcarParser, self).__init__()
        self._filepath = path
        self._filename = filename
        self._parsable_items = OutcarParser.PARSABLE_ITEMS
        self._parsed_data = {}

    def _parse_file(self, inputs):
        """Add all quantities parsed from OUTCAR to _parsed_data."""

        result = self._read_outcar(inputs)
        params = DataFactory('parameter')()
        params.update_dict(result)
        result['parameters'] = params
        return result

    def _read_outcar(self, inputs):
        """Parse the OUTCAR file into a dictionary."""
        result = inputs.get('settings', {})
        result = {}
        energy_free = []
        energy_zero = []
        with open(self._filepath, 'r') as outcar_file_object:
            for line in outcar_file_object:
                # volume
                if line.rfind('volume of cell :') > -1:
                    result['volume'] = float(line.split()[-1])
                # Free energy
                if line.lower().startswith('  free  energy   toten'):
                    energy_free.append(float(line.split()[-2]))
                # Extrapolated zero point energy
                if line.startswith('  energy  without entropy'):
                    energy_zero.append(float(line.split()[-1]))
                # Fermi energy
                if line.rfind('E-fermi') > -1:
                    result['efermi'] = float(line.split()[2])
        result['energies'] = {}
        result['energies']['free_energy'] = energy_free[-1]
        result['energies']['energy_without_entropy'] = energy_zero[-1]
        result['energies']['free_energy_all'] = energy_free
        result['energies']['energy_without_entropy_all'] = energy_zero
        return result
