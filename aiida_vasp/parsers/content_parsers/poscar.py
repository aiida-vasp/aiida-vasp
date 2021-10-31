"""
POSCAR/CONTCAR parser.

--------------
Contains the parsing interfaces to parsevasp used to parse POSCAR/CONTCAR.
"""
# pylint: disable=no-self-use
import numpy as np

from parsevasp.poscar import Poscar, Site
from aiida.common.constants import elements
from aiida_vasp.parsers.content_parsers.parser import BaseFileParser
from aiida_vasp.parsers.node_composer import NodeComposer, get_node_composer_inputs_from_object_parser
from aiida_vasp.utils.aiida_utils import get_data_class

DEFAULT_OPTIONS = {'quantities_to_parse': ['poscar-structure']}


class PoscarParser(BaseFileParser):
    """The parser interface that enables parsing of POSCAR/CONTCAR files.

    The parser is triggered by using the `poscar-structure` quantity key. The quantity key `structure`
    will parse the structure using the XML parser.

    Parameters
    ----------
    precision : int, optional
        An integer specifying the number of digits for floating point numbers that will be written
        to POSCAR/CONTCAR. Defaults to 12.

    """

    PARSABLE_QUANTITIES = {
        'poscar-structure': {
            'inputs': [],
            'name': 'structure',
            'prerequisites': [],
        },
    }

    def __init__(self, *, precision=12, **kwargs):
        """Initialize an instance of this class."""

        self._precision = precision
        super(PoscarParser, self).__init__(**kwargs)

    def _init_from_handler(self, handler):
        """Initialize using a file like handler."""

        try:
            self._content_parser = Poscar(file_handler=handler, prec=self._precision, conserve_order=True, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    def _init_from_data(self, data):
        """Initialize using AiiDA StructureData."""

        if isinstance(data, get_data_class('structure')):
            self._content_data = data
        else:
            raise TypeError('The supplied AiiDA data structure is not a StructureData.')

    def _parse_content(self):
        """Parse the quantities configured and parseable from the content."""

        quantities_to_parse = DEFAULT_OPTIONS.get('quantities_to_parse')
        if self._settings is not None and self._settings.quantity_names_to_parse:
            quantities_to_parse = self._settings.quantity_names_to_parse

        result = {}

        if self._content_parser is None:
            # Parsevasp threw an exception, which means POSCAR could not be parsed.
            for quantity in quantities_to_parse:
                if quantity in self._parsable_quantities:
                    result[quantity] = None
            return result

        for quantity in quantities_to_parse:
            if quantity in self._parsable_quantities:
                # In case there is a - in the quantity, we assume we can
                # parse this quantity from multiple sources, remove source as we do not want to used
                # the source in the property name, i.e. use last element in the split
                quantity_splitted = quantity.split('-')
                quantity_splitted = quantity_splitted[-1]
                result[quantity] = getattr(self, quantity_splitted)

        return result

    @property
    def structure(self):
        return parsevasp_to_aiida(self._content_parser)

    def _content_data_to_content_parser(self):
        """Convert an AiiDA data structure to a content parser instance parsevasp."""
        dictionary = {}
        dictionary['comment'] = self._content_data.label or self._content_data.get_formula()
        dictionary['unitcell'] = np.asarray(self._content_data.cell)
        # As for now all Aiida-structures are in Cartesian coordinates.
        direct = False
        sites = []
        _transform_to_bool = np.vectorize(self.transform_to_bool)
        for index, site in enumerate(self._content_data.sites):
            if self._options is None:
                _selective = [True, True, True]
            else:
                try:
                    _selective = _transform_to_bool(np.array(self._options['positions_dof'])[index, :])
                except KeyError:
                    _selective = [True, True, True]
            sites.append(Site(site.kind_name, site.position, selective=_selective, direct=direct, logger=self._logger))

        dictionary['sites'] = sites

        return Poscar(poscar_dict=dictionary, prec=self._precision, conserve_order=True, logger=self._logger)

    def transform_to_bool(self, value):
        """Helper function to transform the dictionary from strings or integers to bools"""
        if value in [0, 'F', 'f']:
            return False
        if value in [1, 'T', 't']:
            return True
        return True


def parsevasp_to_aiida(poscar):
    """Parsevasp to AiiDA conversion.

    Generate an AiiDA structure from the parsevasp instance of the
    Poscar class.

    Parameters
    ----------
    poscar : object
        An instance of the Poscar class in parsevasp.

    Returns
    -------
    poscar_dict : dict
        A dictionary representation which are ready to be used when creating an
        AiiDA StructureData instance.

    """

    # Fetch a dictionary containing the entries, make sure all coordinates are
    # cartesian
    poscar_dict = poscar.get_dict(direct=False)

    for site in poscar_dict['sites']:
        specie = site['specie']
        # User can specify whatever they want for the elements, but
        # the symbols entries in AiiDA only support the entries defined
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

    return poscar_dict


def fetch_symbols_from_elements(elmnts):
    """Fetch the symbol entry in the elements dictionary in Aiida."""

    new_dict = {}
    for key, value in elmnts.items():
        new_dict[value['symbol']] = key
    return new_dict
