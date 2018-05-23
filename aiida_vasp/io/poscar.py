# pylint: disable=no-self-use
"""Tools for parsing POSCAR files."""
from itertools import groupby

import numpy as np

from parsevasp.poscar import Poscar, Site
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
            return Poscar(poscar_dict=aiida_to_parsevasp(self._data_obj), prec=self.precision, conserve_order=True)
        except SystemExit:
            return None

    def _parse_file(self, inputs):
        """Read POSCAR format file for output structure."""

        result = inputs
        result = {}

        if isinstance(self._data_obj, get_data_class('structure')):
            return {'structure': self._data_obj}

        # parsevasp currently only allows for POSCARs in direct format, while Aiida StructureData
        # only accepts cartesian coordinates. Therefore we have to prepare POSCAR prior to parsing.
        poscar_string, res_dict = prepare_poscar(self._data_obj.path)
        try:
            poscar_dict = Poscar(poscar_string=poscar_string, prec=self.precision, conserve_order=True).entries
            poscar_dict.update(res_dict)
        except SystemExit:
            return {'structure': None}

        result['structure'] = get_data_class('structure')(cell=poscar_dict['unitcell'])

        for site in poscar_dict['sites']:

            position = site.get_position()
            if not res_dict['cartesian']:
                # Convert coordinates to cartesian for Aiida.
                position_cart = poscar_dict['unitcell'][0] * position[0] \
                                + poscar_dict['unitcell'][1] * position[1] \
                                + poscar_dict['unitcell'][2] * position[2]
                position = position_cart

            result['structure'].append_atom(
                position=np.array(position), symbols=res_dict['mapping'][site.get_specie()], name=site.get_specie())

        return result

    def write(self, filepath):
        """
        Override write a POSCAR file from data_obj.

        Since parsevasp modifies the comment in the first line we have to rewrite it one more time.
        """

        if self._parsed_object is not None:
            self._parsed_object.write(filepath)

        with open(filepath, 'r') as file_obj:
            lines = file_obj.readlines()

        # Check whether non-standard types are used and if so replace the comment line.
        lines[0] = self._get_comment_from_lines(lines) + '\n'

        with open(filepath, 'w') as file_obj:
            file_obj.writelines(lines)

    def _get_comment_from_lines(self, lines):
        """
        Prepare the comment line for POSCAR.

        in case non standard kind_names are used, return the mapping kind_name -> symbol.
        Otherwise return the original comment line.
        """

        comment = lines[0].strip()

        # check whether non standard types are used
        mapping = {}
        non_standard_symbol = False

        for site in self._data_obj.sites:
            symbol = self._data_obj.get_kind(site.kind_name).symbols[0]
            mapping[site.kind_name] = site.kind_name
            if site.kind_name != symbol:
                mapping[site.kind_name] = symbol
                non_standard_symbol = True

        if non_standard_symbol:
            # Generate the comment for the mapping of non standard kind_names.
            comment = '# Aiida-elements:'
            for kind_name in lines[5].split():
                comment += ' ' + mapping[kind_name]

        return comment

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
        lines[0] = self._get_comment_from_lines(lines)
        out_string = '\n'.join(lines)

        return out_string


def prepare_poscar(path):
    """
    Prepare POSCAR for parsing with parsevasp.

    If the POSSCAR is in Cartesian format, it has to be converted to Direct coordinates.
    In addition in case that non-standard kind_names have been used, the mapping from
    kind_names to symbols will be read from the comment line.
    """

    result = {}

    with open(path, 'r') as file_obj:
        lines = file_obj.readlines()

    comment = lines[0]
    kind_names = lines[5].split()
    symbols = kind_names

    # Check for non-standard kind_names mapping in the comment line.
    if comment.startswith('# Aiida-elements:'):
        symbols = comment.split(': ')[1].split()

    mapping = {}

    for i, key in enumerate(kind_names):
        mapping[key] = symbols[i]

    result['mapping'] = mapping

    out_string = ''
    cartesian = False

    for line in lines:
        if line.lower().startswith('c'):
            # Change from 'Cartesian' to 'Direct'.
            line = 'Direct\n'
            cartesian = True

        out_string += line

    result['cartesian'] = cartesian

    return out_string, result


def aiida_to_parsevasp(structure):
    """Convert Aiida StructureData to parsevasp's dictionary format."""
    dictionary = {}
    dictionary["comment"] = structure.label or structure.get_formula()
    dictionary["unitcell"] = np.asarray(structure.cell)
    selective = [True, True, True]
    # As for now all Aiida-structures are in Cartesian coordinates.
    direct = True
    sites = []

    for site in structure.sites:

        position = get_direct_coords(dictionary['unitcell'], site.position)
        sites.append(Site(site.kind_name, position, selective=selective, direct=direct))

    dictionary["sites"] = sites
    return dictionary


def get_direct_coords(cell, position):
    """Convert cartesian coordinates to direct coordinates."""

    # First calculate the reciprocal basis.
    omega = np.dot(cell[0], np.cross(cell[1], cell[2]))
    b_0 = 1.0 / omega * np.cross(cell[1], cell[2])
    b_1 = 1.0 / omega * np.cross(cell[2], cell[0])
    b_2 = 1.0 / omega * np.cross(cell[0], cell[1])

    rbasis = np.vstack([b_0, b_1, b_2])

    # Convert to direct coordinates:
    direct = np.asarray([np.dot(position, rbasis[0]), np.dot(position, rbasis[1]), np.dot(position, rbasis[2])])

    return direct
