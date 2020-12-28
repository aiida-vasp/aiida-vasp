"""
POSCAR parser.

--------------
The file parser that handles the parsing of POSCAR and CONTCAR files.
"""
# pylint: disable=no-self-use
import numpy as np

from aiida.common.constants import elements
from parsevasp.poscar import Poscar, Site
from aiida_vasp.parsers.file_parsers.parser import BaseFileParser
from aiida_vasp.parsers.node_composer import NodeComposer
from aiida_vasp.utils.aiida_utils import get_data_class


class PoscarParser(BaseFileParser):
    """
    Parse a POSCAR format file into a StructureData node and vice versa.

    This is a wrapper for parsevasps Poscar class for parsing POSCAR format
    files. The Parsing direction depends on whether the Parser is initialised with
    'path = ...' or 'data = ...'.

    :keyword file_path: Path to the POSCAR file.
    :keyword data: Aiida StructureData for parsing.
    :keyword precision: 'int' number specifying the number of digits for floating point
                        numbers that will be written to POSCAR.
                        DEFAULT = 12

    The PoscarParser will deal with non standard atomic symbols internally if it is
    initialised with StructureData. In case that a POSCAR with non standard atomic
    symbols should be parsed, the comment line must contain the keyword '# Aiida-elements:'
    followed by a list of the actual atomic symbols, e.g.:

        # Aiida-elements: Ga In As
    """

    PARSABLE_ITEMS = {
        'poscar-structure': {
            'inputs': [],
            'name': 'structure',
            'prerequisites': [],
        },
    }

    def __init__(self, *args, **kwargs):
        super(PoscarParser, self).__init__(*args, **kwargs)
        self.precision = 12
        self._structure = None
        self.options = None
        self.init_with_kwargs(**kwargs)

    def _init_with_precision(self, precision):
        self.precision = precision

    def _init_with_options(self, options):
        self.options = options

    def _init_with_data(self, data):
        """Initialize with an AiiDA StructureData instance."""
        if isinstance(data, get_data_class('structure')):
            self._data_obj = data
        else:
            self._logger.warning('Please supply an AiiDA StructureData datatype for `data`.')
            self._data_obj = None
        self._structure = data
        self.parsable_items = self.__class__.PARSABLE_ITEMS
        self._parsed_data = {}

    @property
    def _parsed_object(self):
        """Return the parsevasp object representing the POSCAR file."""

        if isinstance(self._data_obj, get_data_class('structure')):
            # _data_obj is StructureData, return the parsed version if possible.
            try:
                return Poscar(poscar_dict=self.aiida_to_parsevasp(self._data_obj, options=self.options),
                              prec=self.precision,
                              conserve_order=True,
                              logger=self._logger)
            except SystemExit:
                return None
        # _data_obj is a SingleFile:
        return self._data_obj

    def _parse_file(self, inputs):
        """Read POSCAR file format."""

        # check if structure has already been loaded, in that case just return
        if isinstance(self._data_obj, get_data_class('structure')):
            return {'poscar-structure': self._data_obj}

        # pass file path to parsevasp and try to load file
        try:
            poscar = Poscar(file_path=self._data_obj.path, prec=self.precision, conserve_order=True, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally. ' 'Returning None.')
            return {'poscar-structure': None}

        result = parsevasp_to_aiida(poscar)

        return result

    @property
    def structure(self):
        if self._structure is None:
            composer = NodeComposer(file_parsers=[self])
            self._structure = composer.compose('structure', quantities=['poscar-structure'])
        return self._structure

    def aiida_to_parsevasp(self, structure, options=None):
        """Convert Aiida StructureData to parsevasp's dictionary format."""
        dictionary = {}
        dictionary['comment'] = structure.label or structure.get_formula()
        dictionary['unitcell'] = np.asarray(structure.cell)
        # As for now all Aiida-structures are in Cartesian coordinates.
        direct = False
        sites = []
        _transform_to_bool = np.vectorize(self.transform_to_bool)
        for index, site in enumerate(structure.sites):
            if options is None:
                _selective = [True, True, True]
            else:
                try:
                    _selective = _transform_to_bool(np.array(options['positions_dof'])[index, :])
                except KeyError:
                    _selective = [True, True, True]
            sites.append(Site(site.kind_name, site.position, selective=_selective, direct=direct, logger=self._logger))

        dictionary['sites'] = sites
        return dictionary

    def transform_to_bool(self, value):
        """Helper function to transform the dictionary from strings or integers to bools"""
        if value in [0, 'F', 'f']:
            return False
        if value in [1, 'T', 't']:
            return True
        return True


def parsevasp_to_aiida(poscar):
    """
    Parsevasp to Aiida conversion.

    Generate an Aiida structure from the parsevasp instance of the
    Poscar class.

    """

    # fetch a dictionary containing the entries, make sure all coordinates are
    # cartesian
    poscar_dict = poscar.get_dict(direct=False)

    # generate Aiida StructureData and add results from the loaded file
    result = {}

    for site in poscar_dict['sites']:
        specie = site['specie']
        # user can specify whatever they want for the elements, but
        # the symbols entries in Aiida only support the entries defined
        # in aiida.common.constants.elements{}

        # strip trailing _ in case user specifies potential
        symbol = specie.split('_')[0].capitalize()
        # check if leading entry is part of
        # aiida.common.constants.elements{}, otherwise set to X, but first
        # invert
        symbols = fetch_symbols_from_elements(elements)
        try:
            symbols[symbol]
        except KeyError:
            symbol = 'X'

        site['symbol'] = symbol
        site['kind_name'] = specie

    result['poscar-structure'] = poscar_dict

    return result


def fetch_symbols_from_elements(elmnts):
    """Fetch the symbol entry in the elements dictionary in Aiida."""

    new_dict = {}
    for key, value in elmnts.items():
        new_dict[value['symbol']] = key
    return new_dict
