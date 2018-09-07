"""Tools for parsing OUTCAR files."""
import re

from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.io.parser import BaseFileParser

DEFAULT_OPTIONS = {'quantities_to_parse': ['volume', 'energies', 'fermi_level']}


class OutcarParser(BaseFileParser):
    """
    Parse OUTCAR into a dictionary, which is supposed to be turned into ParameterData later.

    For constructor params and more details check the documentation for ``aiida_vasp.io.parser`` and
    ``aiida_vasp.io.parser.BaseParser``.
    """

    FILE_NAME = 'OUTCAR'
    PARSABLE_ITEMS = {
        'outcar-volume': {
            'inputs': ['parameters'],
            'nodeName': 'parameters',
            'prerequisites': []
        },
        'outcar-energies': {
            'inputs': ['parameters'],
            'nodeName': 'parameters',
            'prerequisites': []
        },
        'outcar-fermi_level': {
            'inputs': ['parameters'],
            'nodeName': 'parameters',
            'prerequisites': []
        },
        'outcar-parameters': {
            'inputs': [],
            'nodeName': 'parameters',
            'prerequisites': []
        },
        'symmetries': {
            'inputs': [],
            'nodeName': 'parameters',
            'prerequisites': []
        }
    }

    SPACE_GROUP_OP_PATTERN = re.compile(r'Found\s*(\d+) space group operations')
    POINT_GROUP_OP_PATTERN = re.compile(r'whereof\s*(\d+) operations')
    POINT_SYMMETRY_PATTERN = re.compile(r'point symmetry (.*?)\s*\.')
    SPACE_GROUP_PATTERN = re.compile(r'space group is (.*?)\s*\.')

    def __init__(self, *args, **kwargs):
        super(OutcarParser, self).__init__(*args, **kwargs)
        self.init_with_kwargs(**kwargs)

    def _parse_file(self, inputs):
        """Add all quantities parsed from OUTCAR to _parsed_data."""

        result = self._read_outcar(inputs)
        params = get_data_class('parameter')(dict=result)
        result['outcar-parameters'] = params
        return result

    @staticmethod
    def _parse_line_regex_once(line, regex, res_dict, key, convert=None):
        """
        Parse ``line`` with regular expression ``regex`` optionally converts the result and stores int in ``res_dict[key]``.

        Does not overwrite ``res_dict[key]`` if it already exists and is not None.
        """
        if res_dict.get(key, None) is None:
            regex_result = re.findall(regex, line)
            if not regex_result:
                res_dict[key] = None
            else:
                result = regex_result[0]
                if convert:
                    result = convert(result)
                res_dict[key] = result

    def _read_outcar(self, inputs):
        """Parse the OUTCAR file into a dictionary."""
        result = inputs.get('settings', {})
        result = {}
        energy_free = []
        energy_zero = []
        symmetries = {}
        with open(self._data_obj.path, 'r') as outcar_file_object:
            for line in outcar_file_object:
                # volume
                if line.rfind('volume of cell :') > -1:
                    result['outcar-volume'] = float(line.split()[-1])
                # Free energy
                if line.lower().startswith('  free  energy   toten'):
                    energy_free.append(float(line.split()[-2]))
                # Extrapolated zero point energy
                if line.startswith('  energy  without entropy'):
                    energy_zero.append(float(line.split()[-1]))
                # Fermi energy
                if line.rfind('E-fermi') > -1:
                    result['outcar-efermi'] = float(line.split()[2])
                # space group operations
                self._parse_line_regex_once(line, self.SPACE_GROUP_OP_PATTERN, symmetries, 'num_space_group_operations', int)
                # point group operations
                self._parse_line_regex_once(line, self.POINT_GROUP_OP_PATTERN, symmetries, 'num_point_group_operations', int)
                # point symmetry
                self._parse_line_regex_once(line, self.POINT_SYMMETRY_PATTERN, symmetries, 'point_symmetry')
                # space group
                self._parse_line_regex_once(line, self.SPACE_GROUP_PATTERN, symmetries, 'space_group')
        result['outcar-energies'] = {}
        result['outcar-energies']['free_energy'] = energy_free[-1]
        result['outcar-energies']['energy_without_entropy'] = energy_zero[-1]
        result['outcar-energies']['free_energy_all'] = energy_free
        result['outcar-energies']['energy_without_entropy_all'] = energy_zero
        result['symmetries'] = symmetries
        return result
