"""
Utilities for converting ASE constraints to VASP selective dynamics.

This module provides functions to convert ASE constraint objects (FixAtoms, FixScaled, FixCartesian)
into VASP's selective dynamics format (positions_dof).
"""

from __future__ import annotations

import warnings
from typing import Union

import numpy as np
from aiida import orm
from aiida.common.exceptions import InputValidationError
from aiida.orm.nodes.data.base import to_aiida_type
from ase import Atoms
from ase.constraints import FixAtoms, FixCartesian, FixScaled


def _is_cell_orthogonal(cell: np.ndarray, tol: float = 1e-6) -> bool:
    """
    Check if cell vectors are aligned with Cartesian x, y, z axes.

    :param cell: Cell vectors as (3, 3) array
    :type cell: np.ndarray
    :param tol: Tolerance for considering off-diagonal elements as zero
    :type tol: float
    :returns: True if cell is orthogonal (aligned with Cartesian axes)
    :rtype: bool
    """
    # Check if cell is diagonal (off-diagonal elements should be ~0)
    # Cell format: [[a_x, a_y, a_z], [b_x, b_y, b_z], [c_x, c_y, c_z]]
    # For orthogonal: a_y, a_z, b_x, b_z, c_x, c_y should be ~0
    off_diagonal = [
        cell[0, 1],
        cell[0, 2],  # a_y, a_z
        cell[1, 0],
        cell[1, 2],  # b_x, b_z
        cell[2, 0],
        cell[2, 1],  # c_x, c_y
    ]
    return all(abs(val) < tol for val in off_diagonal)


def atoms_to_positions_dof(atoms: Atoms) -> np.ndarray | None:
    """
    Convert ASE constraints to VASP positions_dof format.

    This function extracts FixAtoms, FixScaled, and FixCartesian constraints from an ASE Atoms object
    and converts them to VASP's selective dynamics format. Multiple constraints are
    properly accumulated using a union of restrictions.

    :param atoms: ASE Atoms object with optional constraints
    :type atoms: Atoms
    :returns: Numpy array of shape (N, 3) with dtype=bool, where True means the atom
        can move in that direction (VASP convention), False means fixed.
        Returns None if no constraints are present.
    :rtype: np.ndarray | None
    :raises InputValidationError: If an unsupported constraint type is encountered or if
        FixCartesian is used with non-orthogonal cells

    .. note::
        - Supports FixAtoms (fixes all 3 directions), FixScaled (fixes specific
          fractional directions), and FixCartesian (fixes Cartesian directions).
        - ASE mask convention (True=fixed) is automatically inverted to VASP
          convention (True=movable).
        - VASP selective dynamics operates in fractional (direct) coordinates.
        - FixCartesian only works correctly when cell vectors are aligned with
          Cartesian axes. A warning is issued if used, and an error is raised if
          the cell is not orthogonal.
        - Multiple constraints on the same atom are accumulated (union of restrictions).

    Example::

        >>> from ase import Atoms
        >>> from ase.constraints import FixAtoms
        >>> atoms = Atoms('H2O', positions=[[0,0,0], [1,0,0], [0,1,0]])
        >>> atoms.set_constraint(FixAtoms(indices=[0]))
        >>> dof = atoms_to_positions_dof(atoms)
        >>> dof[0]  # First atom fixed
        array([False, False, False])
        >>> dof[1]  # Second atom free
        array([True, True, True])
    """
    # Return None if no constraints
    if not atoms.constraints:
        return None

    # Initialize all atoms as movable (True = movable in VASP)
    natoms = len(atoms)
    positions_dof = np.ones((natoms, 3), dtype=bool)

    # Process each constraint
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            # FixAtoms: fix all 3 directions for specified atoms
            indices = constraint.get_indices()
            positions_dof[indices, :] = False

        elif isinstance(constraint, FixScaled):
            # FixScaled: fix specific fractional directions
            # CRITICAL: ASE mask convention is opposite to VASP
            # ASE: True = fixed, False = movable
            # VASP: True = movable, False = fixed
            indices = constraint.get_indices()
            mask = constraint.mask  # ASE convention: True = fixed

            # Invert mask: where ASE mask is True (fixed), set VASP dof to False (fixed)
            positions_dof[indices, mask] = False

        elif isinstance(constraint, FixCartesian):
            # FixCartesian: fix specific Cartesian directions
            # WARNING: VASP selective dynamics works in fractional coordinates!
            # This only works correctly if cell vectors are aligned with Cartesian axes

            # Check if cell is orthogonal
            cell = atoms.get_cell()
            if not _is_cell_orthogonal(cell):
                raise InputValidationError(
                    'FixCartesian constraint cannot be used with non-orthogonal cells. '
                    'VASP selective dynamics operates in fractional coordinates. '
                    'For non-orthogonal cells, use FixScaled instead or ensure your '
                    'cell vectors are aligned with x, y, z axes.'
                )

            # Issue warning about fractional vs Cartesian
            warnings.warn(
                'FixCartesian constraint is being converted to VASP selective dynamics. '
                'Note that VASP selective dynamics operates in fractional (direct) coordinates. '
                'This conversion assumes your cell vectors are aligned with Cartesian x, y, z axes. '
                'Consider using FixScaled for explicit fractional coordinate control.',
                UserWarning,
                stacklevel=2,
            )

            # Treat like FixScaled with same logic
            indices = constraint.get_indices()
            mask = constraint.mask  # ASE convention: True = fixed
            positions_dof[indices, mask] = False

        else:
            # Unsupported constraint type
            constraint_type = type(constraint).__name__
            raise InputValidationError(
                f"Unsupported constraint type '{constraint_type}'. "
                f'Only FixAtoms, FixScaled, and FixCartesian are supported for VASP selective dynamics. '
                f"For complex constraints, manually specify the 'dynamics' port with 'positions_dof'."
            )

    return positions_dof


def serialize_dynamics(atoms: Union[Atoms, dict]) -> orm.Dict | None:
    """
    Serialize ASE Atoms with constraints to dynamics Dict for builder.dynamics port.

    This is a convenience function that converts ASE constraints to an AiiDA Dict
    node ready for direct use with the builder.dynamics port.

    :param atoms: ASE Atoms object with optional constraints
    :type atoms: Atoms
    :returns: orm.Dict with 'positions_dof' key, or None if no constraints present.
        The Dict is ready for direct assignment: ``builder.dynamics = serialize_dynamics(atoms)``
    :rtype: orm.Dict | None

    .. note::
        Handles FixAtoms, FixScaled, and FixCartesian constraints. Multiple constraints are
        properly accumulated (union of restrictions). FixCartesian requires orthogonal cells.

    Example::

        >>> from ase.build import bulk
        >>> from ase.constraints import FixAtoms
        >>> from aiida_vasp.utils import serialize_dynamics
        >>>
        >>> atoms = bulk('Si', 'diamond', a=5.43).repeat((2, 2, 2))
        >>> atoms.set_constraint(FixAtoms(indices=[0, 1, 2, 3]))
        >>>
        >>> # Direct usage with builder
        >>> builder.dynamics = serialize_dynamics(atoms)
    """
    if not isinstance(atoms, Atoms):
        return to_aiida_type(atoms)
    # Get positions_dof array
    dof = atoms_to_positions_dof(atoms)

    # Return None if no constraints
    if dof is None:
        return None

    # Convert to Dict node (tolist() converts numpy array to JSON-serializable list)
    return orm.Dict(dict={'positions_dof': dof.tolist()})
