"""
Generating (symmetry-reduced) k-point grids using spglib
"""

import ase
import numpy as np
from aiida import orm
from aiida.engine import calcfunction
from spglib import get_ir_reciprocal_mesh


def grid_address_to_recip_coord(points, mesh, is_shift=None):
    """
    Convert grid address to fractional coordinates in the reciprocal space
    """
    if is_shift is None:
        shift = np.array((0, 0, 0))
    else:
        shift = np.array([0.5 if shift else 0 for shift in is_shift])
    return (points + shift) / np.asarray(mesh)


def get_ir_kpoints_and_weights(
    atoms: ase.Atoms, mesh, is_time_reversal=True, symprec=1e-5, is_shift=None, symmetry_reduce=True
):
    """
    Return fractional coordinates of irreducible k-points from a given mesh.
    Note: The current implementation does not support using only time-reversal symmetry.

    :param atoms: An ASE atoms object
    :param mesh: A tuple/list for the meshes in each direction or a single number for kpoint distance
    :param is_time_reversal: Whether to use time-reversal symmetry or not.
    :param symprec: Symmetry precision
    :param is_shift: A tuple/list for the shift of the mesh, use [1, 1, 1] for MP Grid.
    :param use_symmetry: Whether to use symmetry or not. If False, the k-points are not reduced at all.

    :return: A tuple of (kpoints, weights).
    """

    # We are actually using a distance rather than a mesh - convert it to a mesh
    if not isinstance(mesh, (list, tuple, np.ndarray)):
        distance = mesh
        the_cell = np.array(atoms.cell)
        reciprocal_cell = 2.0 * np.pi * np.linalg.inv(the_cell).transpose()
        mesh = [max(int(np.ceil(round(np.linalg.norm(b) / distance, 5))), 1) for b in reciprocal_cell]

    spgcell = (atoms.cell, atoms.get_scaled_positions(), atoms.get_atomic_numbers())
    grid_map_table, grid_address = get_ir_reciprocal_mesh(
        mesh, spgcell, is_time_reversal=is_time_reversal, is_shift=is_shift, symprec=symprec, is_dense=False
    )
    if symmetry_reduce:
        unique_point_idx = np.unique(grid_map_table)
        multi = np.array([np.sum(p == grid_map_table) for p in unique_point_idx])  # Compute the multiplicity
        weights = multi / sum(multi)
        coords = grid_address_to_recip_coord(grid_address[unique_point_idx], mesh, is_shift=is_shift)
    else:
        weights = np.ones(len(grid_address)) / len(grid_address)
        coords = grid_address_to_recip_coord(grid_address, mesh, is_shift=is_shift)
    return coords, weights


@calcfunction
def get_ir_kpoints_data(
    structure: orm.StructureData,
    mesh_or_spacing,
    is_time_reversal=True,
    symprec=1e-5,
    is_shift=None,
    symmetry_reduce=True,
):
    """
    Return fractional coordinates of irreducible k-points from a given mesh.
    Note: The current implementation does not support using only time-reversal symmetry.

    :param atoms: An ASE atoms object
    :param mesh: A tuple/list for the meshes in each direction or a single number for kpoint distance
    :param is_time_reversal: Whether to use time-reversal symmetry or not.
    :param symprec: Symmetry precision
    :param is_shift: A tuple/list for the shift of the mesh, use [1, 1, 1] for MP Grid.
    :param use_symmetry: Whether to use symmetry or not. If False, the k-points are not reduced at all.

    :return: A KpointsData object
    """
    if isinstance(is_shift, orm.List):
        is_shift = is_shift.get_list()
    if isinstance(mesh_or_spacing, orm.List):
        mesh_or_spacing = mesh_or_spacing.get_list()
    elif isinstance(mesh_or_spacing, orm.Float):
        mesh_or_spacing = mesh_or_spacing.value
    coords, weights = get_ir_kpoints_and_weights(
        structure.get_ase(),
        mesh_or_spacing,
        is_time_reversal=is_time_reversal.value,
        symprec=symprec.value,
        is_shift=is_shift,
        symmetry_reduce=symmetry_reduce.value,
    )
    kpt = orm.KpointsData()
    kpt.set_cell_from_structure(structure)
    kpt.set_kpoints(coords, weights=weights)
    return kpt
