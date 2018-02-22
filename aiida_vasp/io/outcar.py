"""
Tools for parsing OUTCAR files
"""
from aiida_vasp.io.parser import BaseParser


class OutcarParser(BaseParser):
    """
    parse OUTCAR into a dictionary, which is supposed to be turned into ParameterData later.
    """

    FILE_NAME = 'OUTCAR'
    PARSABLE_ITEMS = {
        'parameters': {
            'volume': None,
            'energies': ['free_energy', 'free_energy_all', 'energy_without_entropy', 'energy_without_entropy_all'],
            'efermi': None,
        },
    }

    def __init__(self, path, filename):
        super(OutcarParser, self).__init__()
        self._filepath = path
        self._filename = filename
        self._parsable_items = OutcarParser.PARSABLE_ITEMS
        self._cached_data = {}

    def _parse_outcar(self):
        """Parse the OUTCAR file into a dictionary."""
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
        result['free_energy'] = energy_free[-1]
        result['energy_without_entropy'] = energy_zero[-1]
        result['free_energy_all'] = energy_free
        result['energy_without_entropy_all'] = energy_zero
        return result

    def _get_volume(self, inputs):
        """Parse the OUTCAR file and return the cell volume"""
        result = inputs
        result.update(self._get_quantity('volume'))
        return {'parameters': result}

    def _get_energies(self, inputs):
        """Parse the OUTCAR file and return the total energies without entropy as well as free energies"""
        result = inputs
        for quantity in self._parsable_items:
            result.update(self._get_quantity(quantity))
        return {'parameters': result}

    def _get_efermi(self, inputs):
        """Return the fermi energy"""
        result = inputs
        result.update(self._get_quantity('efermi'))
        return {'parameters': result}

    def _get_quantity(self, quantity):
        """Return the requested quantity from the _cached_data. If OUTCAR has not been parsed yet, parse it."""
        if not self._cached_data:
            self._cached_data = self._parse_outcar()

        if quantity not in self._cached_data:
            return {}

        return {quantity: self._cached_data[quantity]}
