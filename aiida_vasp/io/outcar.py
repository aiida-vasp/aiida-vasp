"""
Tools for parsing OUTCAR files
"""


class OutcarParser(object):
    """
    parse OUTCAR into a dictionary, which is supposed to be turned into ParameterData later.
    """

    def __init__(self, fname):
        super(OutcarParser, self).__init__()
        self.outcar_file = fname
        self.properties = ['volume', 'free_energy', 'free_energy_all', 'energy_without_entropy', 'energy_without_entropy_all', 'efermi']
        self._cached_data = {}

    @property
    def output_dict(self):
        """Parse the OUTCAR file and return the parsed values wrapped in a dictionary"""
        output_dict = {}
        for property_name in self.properties:
            # _get_quantity returns a dictionary containing the property if parsing has been successful.
            output_dict.update(self._get_quantity(property_name))
        return output_dict

    def _parse_outcar(self):
        """Parse the OUTCAR file into a dictionary."""
        result = {}
        energy_free = []
        energy_zero = []
        with open(self.outcar_file, 'r') as outcar_file_object:
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

    def _get_quantity(self, quantity):
        """Return the requested quantity from the _cached_data. If OUTCAR has not been parsed yet, parse it."""
        if not self._cached_data:
            self._cached_data = self._parse_outcar()

        if quantity not in self._cached_data:
            return {}

        return {quantity: self._cached_data[quantity]}
