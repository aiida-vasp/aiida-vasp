"""
Tools for parsing POSCAR files.

Contains:
- Poscar Parser
- Poscar writer

Features not in `pymatgen.io.vasp.Poscar`:

 * Site ordering is left alone
 * Element grouping allows for different potentials for sites of same species
 * Exposes order of potentials to be used when concatting the POTCAR
"""
from itertools import groupby

import numpy as np
from py import path as py_path  # pylint: disable=no-name-in-module,no-member

from aiida_vasp.io.parser import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class


class PoscarIo(object):
    """
    Write a POSCAR of the form

    <comment>
    <lattice_constant>>
    <lattice vectors>
    <kind counts vector>
    direct
    <sites>

    where:
     * <comment> = label of the structure, or chemical formula if label is ''
     * <lattice_constant> = 1.0
     * <lattice_vector> = 3x3 matrix of floats
     * <kind counts vector> = numbers of sites that should use the same potential consecutively
     * <sites> = <positions> kind_name
     * <positions> = 3-vector of coordinates (in lattice basis)
    """
    POSCAR_TPL = '{comment}\n1.0\n{lattice}\n{kind_counts}\ndirect\n{positions}'
    LATTICE_ROW_TPL = '{:{float_fmt}} {:{float_fmt}} {:{float_fmt}}'
    POS_ROW_TPL = '{:{float_fmt}} {:{float_fmt}} {:{float_fmt}} {label}'

    def __init__(self, structure, precision=None):
        self.structure = structure
        self.float_format = ''
        if precision:
            self.float_format = '.{}'.format(precision)

    def count_kinds(self):
        """
        Count consecutive sites that should use the same potential.

        :return: [(kind_name, num), ... ]
        """
        kind_name_order = [site.kind_name for site in self.structure.sites]
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
        cell = np.array(self.structure.cell)
        if np.linalg.det(cell) < 0:
            cell = cell * -1
        comment = self.structure.label or self.structure.get_formula()
        lattice = '\n'.join([self.LATTICE_ROW_TPL.format(*row, float_fmt=self.float_format) for row in cell])
        kind_counts = ' '.join([str(count[1]) for count in self.count_kinds()])
        positions = '\n'.join([
            self.POS_ROW_TPL.format(*site.position, float_fmt=self.float_format, label=self.structure.get_kind(site.kind_name).symbol)
            for site in self.structure.sites
        ])
        return self.POSCAR_TPL.format(comment=comment, lattice=lattice, kind_counts=kind_counts, positions=positions)

    def write(self, path):
        destination = py_path.local(path)
        destination.write(self.poscar_str())


class PoscarParser(BaseFileParser):
    """Parse a POSCAR format file into a StructureData node."""

    PARSABLE_ITEMS = {
        'structure': {
            'inputs': [],
            'parsers': ['CONTCAR'],
            'nodeName': 'structure',
            'prerequisites': []
        },
    }

    def __init__(self, path, filename, cls):
        super(PoscarParser, self).__init__(cls)
        self._filepath = path
        self._filename = filename
        self._parsable_items = PoscarParser.PARSABLE_ITEMS
        self._parsable_data = {}

    def _parse_file(self, inputs):
        """Read POSCAR format file for output structure."""
        from ase.io import read

        result = inputs
        result = {}
        result['structure'] = get_data_class('structure')()
        cont = self._filepath
        if not cont:
            return {'structure': None}
        result['structure'].set_ase(read(cont, format='vasp'))

        return result
