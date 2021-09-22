"""
POSCAR parser.

--------------
The parser that handles the parsing of POSCAR and CONTCAR.
"""
# pylint: disable=no-self-use
import numpy as np

from parsevasp.poscar import Poscar, Site
from aiida.common.constants import elements
from aiida_vasp.parsers.object_parsers.parser import BaseFileParser
from aiida_vasp.parsers.node_composer import NodeComposer, get_node_composer_inputs_from_object_parser
from aiida_vasp.utils.aiida_utils import get_data_class


class PoscarParser(BaseFileParser):
    """
    Parse a parsevasp representation of POSCAR into a StructureData node and vice versa.

    This is a wrapper for parsevasps Poscar class for parsing POSCARs.
    The Parsing direction depends on whether the Parser is initialised with
    'handler = ...' or 'data = ...'.

    :keyword handler: A handler to a POSCAR.
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
        """
        Initialize POSCAR parser

        data : StructureData

        """

        super(PoscarParser, self).__init__(*args, **kwargs)
        self._structure = None
        self._options = kwargs.get('options', None)
        self._precision = kwargs.get('precision', 12)
        if 'data' in kwargs:
            self._init_structure(kwargs['data'])

    def _init_structure(self, data):
        """Initialize with an AiiDA StructureData instance."""
        if isinstance(data, get_data_class('structure')):
            self._data_obj = data
        else:
            self._logger.warning('Please supply an AiiDA StructureData datatype for `data`.')
            self._data_obj = None
        self._structure = data
        self._parsable_items = self.__class__.PARSABLE_ITEMS
        self._parsed_data = {}

    @property
    def _parsed_object(self):
        """Return the parsevasp object representing the POSCAR."""

        if isinstance(self._data_obj, get_data_class('structure')):
            # _data_obj is StructureData, return the parsed version if possible.
            try:
                return Poscar(poscar_dict=self.aiida_to_parsevasp(self._data_obj, options=self._options),
                              prec=self._precision,
                              conserve_order=True,
                              logger=self._logger)
            except SystemExit:
                return None
        # _data_obj is a SingleFile:
        return self._data_obj

    def _parse_object(self, inputs):
        """Parse the POSCAR object."""

        # Check if structure has already been loaded, in that case just return
        if isinstance(self._data_obj, get_data_class('structure')):
            return {'poscar-structure': self._data_obj}

        # Pass handler to parsevasp and try to initialize the representation
        try:
            poscar = Poscar(file_handler=self._data_obj.handler, prec=self._precision, conserve_order=True, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally. ' 'Returning None.')
            return {'poscar-structure': None}

        result = parsevasp_to_aiida(poscar)

        return result

    @property
    def structure(self):
        if self._structure is None:
            inputs = get_node_composer_inputs_from_object_parser(self, quantity_keys=['poscar-structure'])
            self._structure = NodeComposer.compose('structure', inputs)
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
    Parsevasp to AiiDA conversion.

    Generate an AiiDA structure from the parsevasp instance of the
    Poscar class.
    """

    # Fetch a dictionary containing the entries, make sure all coordinates are
    # cartesian
    poscar_dict = poscar.get_dict(direct=False)

    # Generate AiiDA StructureData and add results from the parsevasp representation
    # of POSCAR.
    result = {}

    for site in poscar_dict['sites']:
        specie = site['specie']
        # User can specify whatever they want for the elements, but
        # the symbols entries in Aiida only support the entries defined
        # in aiida.common.constants.elements{}

        # Strip trailing _ in case user specifies potential
        symbol = specie.split('_')[0].capitalize()
        # Check if leading entry is part of
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
