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
        self.properties = ['volume', 'energies', 'efermi']

    @property
    def output_dict(self):
        """Parse the OUTCAR file and return the parsed values wrapped in a dictionary"""
        output_dict = {}
        for property_name in self.properties:
            try:
                # call the method corresponding to property_name.
                result = getattr(self, '_read_' + property_name)()
                if result is not None:
                    # read_property_name returns a dictionary if parsing has been successful.
                    output_dict.update(result)
            except AttributeError:
                raise NotImplementedError('The OUTCAR parser does not implemnt _read_{}'.format(property_name))
        return output_dict

    def _read_volume(self):
        """Parse the OUTCAR file and return the cell volume"""
        result = {}
        with open(self.outcar_file, 'r') as outcar_file_object:
            for line in outcar_file_object:
                if line.rfind('volume of cell :') > -1:
                    result['volume'] = float(line.split()[-1])
        return result

    def _read_energies(self):
        """Parse the OUTCAR file and return the total energies without entropy as well as free energies"""
        energy_free = []
        energy_zero = []
        with open(self.outcar_file, 'r') as outcar_file_object:
            for line in outcar_file_object:
                # Free energy
                if line.lower().startswith('  free  energy   toten'):
                    energy_free.append(float(line.split()[-2]))
                # Extrapolated zero point energy
                if line.startswith('  energy  without entropy'):
                    energy_zero.append(float(line.split()[-1]))
        result = {}
        result['free_energy'] = energy_free[-1]
        result['energy_without_entropy'] = energy_zero[-1]
        result['free_energy_all'] = energy_free
        result['energy_without_entropy_all'] = energy_zero
        return result

    def _read_efermi(self):
        """Parse the OUTCAR file and return the fermi energy"""
        result = {}
        with open(self.outcar_file, 'r') as outcar_file_object:
            for line in outcar_file_object:
                if line.rfind('E-fermi') > -1:
                    result['efermi'] = float(line.split()[2])
        return result
