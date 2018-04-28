# pylint: disable=no-self-use
"""Tools for parsing POSCAR files."""
from itertools import groupby

import numpy as np
from py import path as py_path  # pylint: disable=no-name-in-module,no-member

from parsevasp.poscar import Poscar, Site
from aiida_vasp.io.parser import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class


class PoscarParser(BaseFileParser):
    """
    Parse a POSCAR format file into a StructureData node and vice versa.

    This is a wrapper for parsevasps Poscar class for parsing POSCAR format
    files. The Parsing direction depends on whether the Parser is initialised with
    'path = ...' or 'data = ...'.

    The PoscarParser will deal with non standard atomic symbols internally if it is
    initialised with StructureData. In case that a POSCAR with non standard atomic
    symbols should be parsed, the comment line must contain the keyword '# Aiida-elements:'
    and a list of the actual atomic symbols, e.g.:

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

    POSCAR_TPL = '{comment}\n1.0\n{lattice}\n{kind_counts}\ncartesian\n{positions}'
    LATTICE_ROW_TPL = '{:{float_fmt}} {:{float_fmt}} {:{float_fmt}}'
    POS_ROW_TPL = '{:{float_fmt}} {:{float_fmt}} {:{float_fmt}} {label}'

    def __init__(self, *args, **kwargs):
        super(PoscarParser, self).__init__(*args, **kwargs)
        self.float_format = ''
        self.init_with_kwargs(**kwargs)

    def _init_with_precision(self, precision):
        self.float_format = '.{}'.format(precision)

    def _init_with_data(self, data):
        """Init with Aiida StructureData"""
        self._data_obj = data
        self._parsable_items = self.__class__.PARSABLE_ITEMS
        self._parsed_data = {}

    @property
    def _parsed_object(self):
        """Return the parsevasp object representing the POSCAR file."""

        try:
            return Poscar(poscar_dict=aiida_to_parsevasp(self._data_obj))
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
            poscar_dict = Poscar(poscar_string=poscar_string).entries
            poscar_dict.update(res_dict)
        except SystemExit:
            return {'structure': None}

        result['structure'] = get_data_class('structure')(cell=poscar_dict['unitcell'])

        print res_dict['mapping']

        for site in poscar_dict['sites']:

            position = site.get_position()
            if not res_dict['cartesian']:
                position_cart = poscar_dict['unitcell'][0] * position[0] \
                                + poscar_dict['unitcell'][1] * position[1] \
                                + poscar_dict['unitcell'][2] * position[2]
                position = position_cart

            result['structure'].append_atom(
                position=np.array(position), symbols=res_dict['mapping'][site.get_specie()], name=site.get_specie())

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
        return [kind[0] for kind in self.count_kinds()]

    def poscar_str(self):
        """
        Create a string of the POSCAR contents.

        Accounts for lattices which have triple product < 0 by inverting lattice vectors in that case.
        """
        cell = np.array(self._data_obj.cell)
        if np.linalg.det(cell) < 0:
            cell = cell * -1
        comment = self._data_obj.label or self._data_obj.get_formula()
        lattice = '\n'.join([self.LATTICE_ROW_TPL.format(*row, float_fmt=self.float_format) for row in cell])
        kind_counts = ' '.join([str(count[1]) for count in self.count_kinds()])
        positions = '\n'.join([
            self.POS_ROW_TPL.format(*site.position, float_fmt=self.float_format, label=self._data_obj.get_kind(site.kind_name).symbol)
            for site in self._data_obj.sites
        ])
        return self.POSCAR_TPL.format(comment=comment, lattice=lattice, kind_counts=kind_counts, positions=positions)

    def write(self, filepath):
        destination = py_path.local(filepath)
        destination.write(self.poscar_str())


def prepare_poscar(path):
    """Workaround for Cartesian coordinate POSCAR, which is currently not supported by parsevasp."""

    result = {}

    with open(path, 'r') as file_obj:
        lines = file_obj.readlines()

    comment = lines[0]
    kind_names = lines[5].split()
    symbols = kind_names

    print kind_names

    if comment.startswith('# Aiida-elements:'):
        symbols = comment.split(':')[1].split(' ')

    mapping = {}

    for i, key in enumerate(kind_names):
        mapping[key] = symbols[i]

    result['mapping'] = mapping

    out_string = ""
    cartesian = False
    for line in lines:
        if line.lower().startswith('c'):
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
        # position = get_direct_coords(dictionary['unitcell'], np.asarray(site.position))
        position = np.asarray(site.position)
        sites.append(Site(site.kind_name, position, selective=selective, direct=direct))
    dictionary["sites"] = sites
    return dictionary


def get_direct_coords(cell, position):
    """Convert cartesian coordinates to direct coordinates."""

    # First calculate the reciprocal basis

    omega = np.dot(cell[0], np.cross(cell[1], cell[2]))
    b_0 = 1.0 / omega * np.cross(cell[1], cell[2])
    b_1 = 1.0 / omega * np.cross(cell[2], cell[0])
    b_2 = 1.0 / omega * np.cross(cell[0], cell[1])

    rbasis = np.asarray([b_0, b_1, b_2])

    # Convert to direct coordinates:

    direct = np.asarray([np.dot(position, rbasis[0]), np.dot(position, rbasis[1])], np.dot(position, rbasis[2]))

    return direct
