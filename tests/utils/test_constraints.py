"""Tests for constraints module."""

import numpy as np
import pytest
from aiida import orm
from aiida.common.exceptions import InputValidationError
from ase import Atoms
from ase.constraints import FixAtoms, FixBondLength, FixCartesian, FixScaled

from aiida_vasp.utils.constraints import atoms_to_positions_dof, serialize_dynamics
from aiida_vasp.workchains.v2.relax import get_maximum_force


def test_fixatoms_basic():
    """Test FixAtoms constraint fixes all 3 directions for specified atoms."""
    atoms = Atoms('H5', positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0], [4, 0, 0]])
    atoms.set_constraint(FixAtoms(indices=[0, 2, 4]))

    dof = atoms_to_positions_dof(atoms)

    expected = np.array(
        [
            [False, False, False],  # atom 0: fixed
            [True, True, True],  # atom 1: free
            [False, False, False],  # atom 2: fixed
            [True, True, True],  # atom 3: free
            [False, False, False],  # atom 4: fixed
        ]
    )

    np.testing.assert_array_equal(dof, expected)


def test_fixscaled_basic():
    """Test FixScaled with partial constraints."""
    atoms = Atoms('H3', positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]], cell=[5, 5, 5])
    # Fix x and z directions for atom 1 (ASE: True=fixed)
    atoms.set_constraint(FixScaled([1], mask=(True, False, True)))

    dof = atoms_to_positions_dof(atoms)

    expected = np.array(
        [
            [True, True, True],  # atom 0: free
            [False, True, False],  # atom 1: x and z fixed, y movable
            [True, True, True],  # atom 2: free
        ]
    )

    np.testing.assert_array_equal(dof, expected)


def test_fixscaled_mask_inversion():
    """Test that ASE mask convention (True=fixed) is inverted to VASP (True=movable).

    This is the CRITICAL test for mask inversion.
    """
    atoms = Atoms('H2', positions=[[0, 0, 0], [1, 0, 0]], cell=[5, 5, 5])

    # ASE convention: mask=(True, True, False) means fix x and y, allow z
    atoms.set_constraint(FixScaled([1], mask=(True, True, False)))

    dof = atoms_to_positions_dof(atoms)

    # VASP convention: (False, False, True) means x and y fixed, z movable
    expected = np.array(
        [
            [True, True, True],  # atom 0: all movable
            [False, False, True],  # atom 1: x,y fixed (inverted from ASE), z movable
        ]
    )

    np.testing.assert_array_equal(dof, expected)


def test_multiple_fixatoms():
    """Test multiple FixAtoms constraints are properly accumulated."""
    atoms = Atoms('H4', positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])

    # Multiple FixAtoms constraints
    constraint1 = FixAtoms(indices=[0])
    constraint2 = FixAtoms(indices=[2, 3])
    atoms.set_constraint([constraint1, constraint2])

    dof = atoms_to_positions_dof(atoms)

    expected = np.array(
        [
            [False, False, False],  # atom 0: fixed
            [True, True, True],  # atom 1: free
            [False, False, False],  # atom 2: fixed
            [False, False, False],  # atom 3: fixed
        ]
    )

    np.testing.assert_array_equal(dof, expected)


def test_multiple_mixed():
    """Test FixAtoms + FixScaled constraints together."""
    atoms = Atoms('H4', positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], cell=[5, 5, 5])

    # Fix atom 0 completely, fix z-direction for atoms 2 and 3
    constraint1 = FixAtoms(indices=[0])
    constraint2 = FixScaled([2, 3], mask=(False, False, True))
    atoms.set_constraint([constraint1, constraint2])

    dof = atoms_to_positions_dof(atoms)

    expected = np.array(
        [
            [False, False, False],  # atom 0: fully fixed (FixAtoms)
            [True, True, True],  # atom 1: free
            [True, True, False],  # atom 2: z fixed (FixScaled)
            [True, True, False],  # atom 3: z fixed (FixScaled)
        ]
    )

    np.testing.assert_array_equal(dof, expected)


def test_overlapping_constraints():
    """Test overlapping constraints on same atom (union of restrictions)."""
    atoms = Atoms('H2', positions=[[0, 0, 0], [1, 0, 0]], cell=[5, 5, 5])

    # Fix x direction
    constraint1 = FixScaled([0], mask=(True, False, False))
    # Fix y direction on same atom
    constraint2 = FixScaled([0], mask=(False, True, False))
    atoms.set_constraint([constraint1, constraint2])

    dof = atoms_to_positions_dof(atoms)

    # Union: both x and y should be fixed
    expected = np.array(
        [
            [False, False, True],  # atom 0: x and y fixed, z movable
            [True, True, True],  # atom 1: free
        ]
    )

    np.testing.assert_array_equal(dof, expected)


def test_no_constraints():
    """Test that None is returned when no constraints are present."""
    atoms = Atoms('H3', positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])

    dof = atoms_to_positions_dof(atoms)

    assert dof is None


def test_unsupported_constraint():
    """Test that unsupported constraint types raise InputValidationError."""
    atoms = Atoms('H2', positions=[[0, 0, 0], [1, 0, 0]])
    atoms.set_constraint(FixBondLength(0, 1))

    with pytest.raises(InputValidationError, match='Unsupported constraint type'):
        atoms_to_positions_dof(atoms)


def test_serialize_dynamics(aiida_profile):
    """Test serialize_dynamics returns correct orm.Dict."""
    atoms = Atoms('H3', positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
    atoms.set_constraint(FixAtoms(indices=[0, 2]))

    result = serialize_dynamics(atoms)

    # Check return type
    assert isinstance(result, orm.Dict)

    # Check content
    dof_dict = result.get_dict()
    assert 'positions_dof' in dof_dict

    dof = dof_dict['positions_dof']
    expected = [
        [False, False, False],  # atom 0: fixed
        [True, True, True],  # atom 1: free
        [False, False, False],  # atom 2: fixed
    ]

    assert dof == expected


def test_serialize_no_constraints(aiida_profile):
    """Test serialize_dynamics returns None for no constraints."""
    atoms = Atoms('H3', positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])

    result = serialize_dynamics(atoms)

    assert result is None


def test_serialize_function(aiida_profile):
    """Test serialize_dynamics returns input dicts unchanged."""

    assert serialize_dynamics({'x': 1}).get_dict() == {'x': 1}


def test_empty_atoms():
    """Test graceful handling of empty atoms object."""
    atoms = Atoms()

    dof = atoms_to_positions_dof(atoms)

    assert dof is None


def test_fixscaled_single_direction():
    """Test FixScaled fixing only single direction."""
    atoms = Atoms('H2', positions=[[0, 0, 0], [1, 0, 0]], cell=[5, 5, 5])

    # Fix only z-direction (most common for slabs)
    atoms.set_constraint(FixScaled([0], mask=(False, False, True)))

    dof = atoms_to_positions_dof(atoms)

    expected = np.array(
        [
            [True, True, False],  # atom 0: only z fixed
            [True, True, True],  # atom 1: free
        ]
    )

    np.testing.assert_array_equal(dof, expected)


def test_fixatoms_with_multiple_atoms():
    """Test FixAtoms with many atoms (realistic case)."""
    # Create a larger structure
    atoms = Atoms('H10', positions=[[i, 0, 0] for i in range(10)])

    # Fix bottom 3 atoms (like a slab)
    atoms.set_constraint(FixAtoms(indices=[0, 1, 2]))

    dof = atoms_to_positions_dof(atoms)

    # First 3 should be fixed
    assert np.all(~dof[:3])
    # Rest should be free
    assert np.all(dof[3:])


def test_fixcartesian_orthogonal():
    """Test FixCartesian with orthogonal cell (should work with warning)."""
    # Create orthogonal cell (aligned with x, y, z)
    atoms = Atoms('H3', positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]], cell=[5, 5, 5], pbc=True)

    # Fix z-direction in Cartesian coordinates
    atoms.set_constraint(FixCartesian([1], mask=(False, False, True)))

    # Should work but issue a warning
    with pytest.warns(UserWarning, match='FixCartesian constraint is being converted'):
        dof = atoms_to_positions_dof(atoms)

    expected = np.array(
        [
            [True, True, True],  # atom 0: free
            [True, True, False],  # atom 1: z fixed
            [True, True, True],  # atom 2: free
        ]
    )

    np.testing.assert_array_equal(dof, expected)


def test_fixcartesian_non_orthogonal():
    """Test FixCartesian with non-orthogonal cell (should raise error)."""
    # Create non-orthogonal cell
    atoms = Atoms(
        'H2',
        positions=[[0, 0, 0], [1, 0, 0]],
        cell=[[5, 0, 0], [2, 5, 0], [0, 0, 5]],  # Non-orthogonal: b_x != 0
        pbc=True,
    )

    atoms.set_constraint(FixCartesian([0], mask=(True, False, False)))

    with pytest.raises(InputValidationError, match='FixCartesian constraint cannot be used with non-orthogonal cells'):
        atoms_to_positions_dof(atoms)


def test_fixcartesian_multiple_directions():
    """Test FixCartesian fixing multiple directions."""
    atoms = Atoms('H2', positions=[[0, 0, 0], [1, 0, 0]], cell=[5, 5, 5], pbc=True)

    # Fix x and y directions
    atoms.set_constraint(FixCartesian([0], mask=(True, True, False)))

    with pytest.warns(UserWarning):
        dof = atoms_to_positions_dof(atoms)

    expected = np.array(
        [
            [False, False, True],  # atom 0: x,y fixed, z free
            [True, True, True],  # atom 1: free
        ]
    )

    np.testing.assert_array_equal(dof, expected)


def test_mixed_constraints_with_fixcartesian():
    """Test mixing FixAtoms, FixScaled, and FixCartesian."""
    atoms = Atoms('H5', positions=[[i, 0, 0] for i in range(5)], cell=[10, 10, 10], pbc=True)

    constraint1 = FixAtoms(indices=[0])  # Fully fix atom 0
    constraint2 = FixScaled([1], mask=(False, False, True))  # Fix z for atom 1
    constraint3 = FixCartesian([2], mask=(True, False, False))  # Fix x for atom 2

    atoms.set_constraint([constraint1, constraint2, constraint3])

    with pytest.warns(UserWarning):
        dof = atoms_to_positions_dof(atoms)

    expected = np.array(
        [
            [False, False, False],  # atom 0: fully fixed (FixAtoms)
            [True, True, False],  # atom 1: z fixed (FixScaled)
            [False, True, True],  # atom 2: x fixed (FixCartesian)
            [True, True, True],  # atom 3: free
            [True, True, True],  # atom 4: free
        ]
    )

    np.testing.assert_array_equal(dof, expected)


def test_integration_with_get_maximum_force():
    """Test integration with get_maximum_force from relax workchain.

    This verifies the complete workflow: ASE constraints → positions_dof → force calculation.
    Forces on fixed atoms should be ignored when calculating maximum force.
    """
    # Create atoms with mixed constraints
    atoms = Atoms('H5', positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0], [4, 0, 0]])
    # Fix atoms 0 and 2 completely, fix z-direction for atom 4
    constraint1 = FixAtoms(indices=[0, 2])
    constraint2 = FixScaled([4], mask=(False, False, True))
    atoms.set_constraint([constraint1, constraint2])

    # Convert constraints to positions_dof
    dof = atoms_to_positions_dof(atoms)

    # Expected dof array
    expected_dof = np.array(
        [
            [False, False, False],  # atom 0: fully fixed
            [True, True, True],  # atom 1: free
            [False, False, False],  # atom 2: fully fixed
            [True, True, True],  # atom 3: free
            [True, True, False],  # atom 4: z fixed
        ]
    )
    np.testing.assert_array_equal(dof, expected_dof)

    # Create force array where fixed atoms have large forces
    forces = np.array(
        [
            [10.0, 10.0, 10.0],  # atom 0: large force but FIXED → should be ignored
            [0.1, 0.2, 0.3],  # atom 1: small force, free
            [5.0, 5.0, 5.0],  # atom 2: large force but FIXED → should be ignored
            [0.4, 0.5, 0.6],  # atom 3: medium force, free
            [1.0, 1.0, 1.0],  # atom 4: medium force but z is FIXED → should be ignored
        ]
    )

    # Without dof: maximum force would be from atom 0 (norm = 10*sqrt(3) ≈ 17.32)
    max_force_no_dof = get_maximum_force(forces)
    assert max_force_no_dof == pytest.approx(10.0 * np.sqrt(3), rel=1e-6)

    # With dof: atoms 0, 2, and 4 are ignored, maximum is from atom 3
    max_force_with_dof = get_maximum_force(forces, dof=dof)
    expected_max = np.sqrt(0.4**2 + 0.5**2 + 0.6**2)
    assert max_force_with_dof == pytest.approx(expected_max, rel=1e-6)

    # Verify forces on fixed atoms are properly ignored
    assert max_force_with_dof < 1.0  # Much smaller than forces on fixed atoms
