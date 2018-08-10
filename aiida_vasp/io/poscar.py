# pylint: disable=no-self-use
"""Tools for parsing POSCAR files."""
import numpy as np

from parsevasp.poscar import Poscar, Site
from aiida.common.constants import elements
from aiida_vasp.io.parser import BaseFileParser
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
            'nodeName': 'structure',
            'prerequisites': [],
            'alternatives': ['structure']
        },
    }

    def __init__(self, *args, **kwargs):
        super(PoscarParser, self).__init__(*args, **kwargs)
        self.precision = 12
        self.init_with_kwargs(**kwargs)

    def _init_with_precision(self, precision):
        self.precision = precision

    def _init_with_data(self, data):
        """Init with Aiida StructureData"""
        self._data_obj = data
        self._parsable_items = self.__class__.PARSABLE_ITEMS
        self._parsed_data = {}

    @property
    def _parsed_object(self):
        """Return the parsevasp object representing the POSCAR file."""

        try:
            return Poscar(poscar_dict=aiida_to_parsevasp(self._parsed_data['poscar-structure']), prec=self.precision, conserve_order=True)
        except SystemExit:
            return None

    def _parse_file(self, inputs):
        """Read POSCAR file format."""

        # check if structure have already been loaded, in that case just return
        if isinstance(self._data_obj, get_data_class('structure')):
            return {'poscar-structure': self._data_obj}

        # pass file path to parsevasp and try to load file
        try:
            poscar = Poscar(file_path=self._data_obj.path, prec=self.precision, conserve_order=True)
        except SystemExit:
            self._logger.warning("Parsevasp exited abnormally. " "Returning None.")
            return {'poscar-structure': None}

        result = parsevasp_to_aiida(poscar)

        return result


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

    result['poscar-structure'] = get_data_class('structure')(cell=poscar_dict['unitcell'])

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
        result['poscar-structure'].append_atom(position=site['position'], symbols=symbol, name=specie)

    return result


def aiida_to_parsevasp(structure):
    """Convert Aiida StructureData to parsevasp's dictionary format."""
    dictionary = {}
    dictionary["comment"] = structure.label or structure.get_formula()
    dictionary["unitcell"] = np.asarray(structure.cell)
    selective = [True, True, True]
    # As for now all Aiida-structures are in Cartesian coordinates.
    direct = False
    sites = []
    for site in structure.sites:
        sites.append(Site(site.kind_name, site.position, selective=selective, direct=direct))

    dictionary["sites"] = sites
    return dictionary


def fetch_symbols_from_elements(elmnts):
    """Fetch the symbol entry in the elements dictionary in Aiida."""

    new_dict = {}
    for key, value in elmnts.items():
        new_dict[value['symbol']] = key
    return new_dict
