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
        self.init_with_kwargs(**kwargs)

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

        poscar_string, cartesian = prepare_poscar(self._data_obj.path)
        try:
            poscar_dict = Poscar(poscar_string=poscar_string).entries
            poscar_dict['cartesian'] = cartesian
        except SystemExit:
            return {'structure': None}

        result['structure'] = get_data_class('structure')(cell=poscar_dict['unitcell'])
        # parsevasp currently only allows for POSCARs in direct format, while Aiida StructureData
        # only accepts cartesian coordinates.
        for site in poscar_dict['sites']:

            position = site.get_position()
            if not cartesian:
                position_cart = poscar_dict['unitcell'][0] * position[0] \
                                + poscar_dict['unitcell'][1] * position[1] \
                                + poscar_dict['unitcell'][2] * position[2]
                position = position_cart

            result['structure'].append_atom(position=np.array(position), symbols=site.get_specie())

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


def prepare_poscar(path):
    """Workaround for Cartesian coordinate POSCAR, which is currently not supported by parsevasp."""

    with open(path, 'r') as file_obj:
        lines = file_obj.readlines()

    out_string = ""
    cartesian = False
    for line in lines:
        if line.lower().startswith('c'):
            line = 'Direct'
            cartesian = True
        out_string += line

    return out_string, cartesian


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
        position = get_direct_coords(dictionary['unitcell'], np.asarray(site.position))
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
