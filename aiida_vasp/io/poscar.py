# pylint: disable=no-self-use
"""Tools for parsing POSCAR files."""

import numpy as np

from parsevasp.poscar import Poscar, Site
from aiida_vasp.io.parser import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class


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
        position = np.asarray(site.position)
        sites.append(Site(site.kind_name, position, selective=selective, direct=direct))
    dictionary["sites"] = sites
    return dictionary


class PoscarParser(BaseFileParser):
    """
    Parse a POSCAR format file into a StructureData node and vice versa.

    The Parsing direction depends on whether the Parser is initialised with
    'path = ...' or 'structure = ...'.
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

    def _init_with_path(self, path):
        self._filepath = path
        self._parsable_items = PoscarParser.PARSABLE_ITEMS
        self._parsed_data = {}

    def _init_with_structure(self, structure):
        """Init with Aiida StructureData"""
        self._data_obj = structure

        try:
            self._parsed_obj = Poscar(poscar_dict=aiida_to_parsevasp(structure))
        except SystemExit:
            self._parsed_obj = None

    def _parse_file(self, inputs):
        """Read POSCAR format file for output structure."""

        result = inputs
        result = {}

        try:
            poscar_dict = Poscar(file_path=self._filepath).entries
        except SystemExit:
            return {'structure': None}

        result['structure'] = get_data_class('structure')(cell=poscar_dict['unitcell'])
        # parsevasp currently only allows for POSCARs in direct format, while Aiida StructureData
        # only accepts cartesian coordinates.
        for site in poscar_dict['sites']:
            result['structure'].append(position=np.array(site.get_position()), symbols=site.get_specie())

        return result
