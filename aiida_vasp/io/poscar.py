# pylint: disable=no-self-use
"""Tools for parsing POSCAR files."""
from itertools import groupby
from collections import Counter, defaultdict
import numpy as np
import sys


from parsevasp.poscar import Poscar, Site
from aiida_vasp.io.parser import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida.common.constants import elements

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
        'structure': {
            'inputs': [],
            'parsers': ['CONTCAR'],
            'nodeName': 'structure',
            'prerequisites': []
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
            return Poscar(poscar_dict=aiida_to_parsevasp(self._parsed_data['structure']),
                          prec=self.precision, conserve_order=True)
        except SystemExit:
            return None

    def _parse_file(self, inputs):
        """Read POSCAR file format."""

        # check if structure have already been loaded, in that case just return
        if isinstance(self._data_obj, get_data_class('structure')):
            return {'structure': self._data_obj}

        # pass file path to parsevasp and try to load file
        try:
            poscar = Poscar(file_path=self._data_obj.path, prec=self.precision,
                            conserve_order=True)
        except SystemExit:
            self._logger.warning("Parsevasp exited abnormally. "
                                 "Returning None.")
            return {'structure': None}

        result = parsevasp_to_aiida(poscar)

        return result

    def count_kinds(self):
        """
        Count consecutive sites that should use the same potential.

        :return: [(kind_name, num), ... ]
        """
        kind_name_order = [site.kind_name for site in self._data_obj.sites]
        groups = groupby(kind_name_order)
        counts = [(label, sum(1 for _ in group)) for label, group in groups]
        return counts

    @property
    def potentials_order(self):
        """Get the order of the potentials from the StructureData stored in _data_obj."""
        return [kind[0] for kind in self.count_kinds()]

    def get_string(self):
        """Create a string of the POSCAR contents."""

        lines = self._parsed_object.get_string().split('\n')
        # Parsevasp replaces the comment line, so we have to update it in case
        # that non standard kind_names have been used.
        # eFL : commented this functionality for the moment
        # lines[0] = self._get_comment_from_lines(lines)
        out_string = '\n'.join(lines)

        return out_string

def parsevasp_to_aiida(poscar):
    """Generate an Aiida structure from the parsevasp instance of the
    Poscar class.

    """
    
    # fetch a dictionary containing the entries, make sure all coordinates are
    # cartesian
    poscar_dict = poscar.get_dict(direct = False)
    
    # generate Aiida StructureData and add results from the loaded file
    result = {}
    
    result['structure'] = get_data_class('structure') \
                          (cell=poscar_dict['unitcell'])
    
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
        result['structure'].append_atom(position=site['position'],
                                        symbols=symbol, name=specie)

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
        sites.append(Site(site.kind_name, site.position,
                          selective=selective, direct=direct))

    dictionary["sites"] = sites
    return dictionary

def fetch_symbols_from_elements(elmnts):
    """Fetch the symbol entry in the elements
    dictionary in Aiida.

    """

    new_dict = {}
    for k, v in elmnts.items():
        new_dict[v['symbol']]=k
    return new_dict
