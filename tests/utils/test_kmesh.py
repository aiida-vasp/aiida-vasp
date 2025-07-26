import ase
import numpy as np
import pytest
from aiida import orm

from aiida_vasp.utils.kmesh import get_ir_kpoints_and_weights, get_ir_kpoints_data


@pytest.fixture
def atoms():
    """A test structure"""
    lattice = np.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]) * 5.4
    positions = [[0.875, 0.875, 0.875], [0.125, 0.125, 0.125]]
    numbers = [
        1,
    ] * 2
    return ase.Atoms(cell=lattice, scaled_positions=positions, numbers=numbers)


@pytest.fixture
def structure(atoms):
    return orm.StructureData(ase=atoms)


def test_get_ir_kpoints_from_structure(atoms):
    """Test the get_ir_kpoints_from_structure function"""

    mesh = [8, 8, 8]

    # Test gamma-center grids
    ir_kpoints, weights = get_ir_kpoints_and_weights(atoms.cell, atoms.get_scaled_positions(), atoms.numbers, mesh)
    assert ir_kpoints.shape == (29, 3)
    assert np.sum(weights) == 1
    # Test shifted grids
    ir_kpoints, weights = get_ir_kpoints_and_weights(
        atoms.cell, atoms.get_scaled_positions(), atoms.numbers, mesh, is_shift=[1, 1, 1]
    )
    assert ir_kpoints.shape == (60, 3)
    assert np.sum(weights) == 1

    # Test not doing reduction at all
    ir_kpoints, weights = get_ir_kpoints_and_weights(
        atoms.cell, atoms.get_scaled_positions(), atoms.numbers, mesh, symmetry_reduce=False
    )
    assert ir_kpoints.shape == (8 * 8 * 8, 3)
    assert np.sum(weights) == 1


def test_get_ir_kpoints_data(structure):
    """Test the get_ir_kpoints_data function"""
    mesh = [8, 8, 8]
    kpt = get_ir_kpoints_data(structure, mesh)
    coords, weights = kpt.get_kpoints(also_weights=True)
    assert coords.shape == (29, 3)
    assert sum(weights) == 1

    kpt = get_ir_kpoints_data(structure, mesh, is_shift=[1, 1, 1])
    coords, weights = kpt.get_kpoints(also_weights=True)
    assert coords.shape == (60, 3)
    assert sum(weights) == 1

    kpt = get_ir_kpoints_data(structure, mesh, symmetry_reduce=False)
    coords, weights = kpt.get_kpoints(also_weights=True)
    assert coords.shape == (8 * 8 * 8, 3)
    assert sum(weights) == 1
