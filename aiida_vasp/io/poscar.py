from itertools import groupby

import numpy as np
from py import path as py_path  # pylint: disable=no-name-in-module,no-member


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

    def poscar_str(self):
        """
        Create a string of the POSCAR contents.

        Accounts for lattices which have triple product < 0 by inverting lattice vectors in that case.
        """
        cell = np.array(self.structure.cell)
        if np.linalg.det(cell) < 0:
            cell = -cell
        comment = self.structure.label or self.structure.get_formula()
        lattice = '\n'.join([self.LATTICE_ROW_TPL.format(*row, float_fmt=self.float_format) for row in cell])
        kind_counts = ' '.join([str(count[1]) for count in self.count_kinds()])
        positions = '\n'.join([self.POS_ROW_TPL.format(*site.position, float_fmt=self.float_format, label=site.kind_name) for site in self.structure.sites])
        return self.POSCAR_TPL.format(comment=comment, lattice=lattice, kind_counts=kind_counts, positions=positions)

    def write(self, path):
        destination = py_path.local(path)
        destination.write(self.poscar_str())

