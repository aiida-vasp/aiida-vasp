"""
Tests for :mod:`aiida_vasp.common.magmapping`.
"""

from aiida_vasp.common.magmapping import convert_to_plain_list, create_additional_species


def test_create_additional_species_scalar():
    """Anti-ferromagnetic two sublattice Fe."""
    species = ['Fe', 'O', 'Fe', 'O']
    magmom = [2.0, 0.0, -2.0, 0.0]
    new_species, mapping = create_additional_species(species, magmom)
    # Fe atoms with opposite spins must end up as different species (Fe1, Fe2)
    assert new_species == ['Fe1', 'O', 'Fe2', 'O']
    assert mapping == {'Fe1': 2.0, 'Fe2': -2.0, 'O': 0.0}


def test_create_additional_species_unique_magmom_kept():
    """Single-species magmom should keep the original name."""
    species = ['Fe', 'Fe']
    magmom = [2.0, 2.0]
    new_species, mapping = create_additional_species(species, magmom)
    assert new_species == ['Fe', 'Fe']
    assert mapping == {'Fe': 2.0}


def test_create_additional_species_vector():
    """3D vector magmom should produce distinct species per unique vector."""
    species = ['Fe', 'O', 'Fe', 'O']
    magmom = [(2.0, 0.0, 0.0), (0.0, 0.0, 0.0), (-2.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
    new_species, mapping = create_additional_species(species, magmom)
    assert new_species == ['Fe1', 'O', 'Fe2', 'O']
    assert mapping == {'Fe1': (2.0, 0.0, 0.0), 'Fe2': (-2.0, 0.0, 0.0), 'O': (0.0, 0.0, 0.0)}


def test_create_additional_species_mixed_scalar_vector():
    """Mixed scalar and vector magmoms must be treated as different species buckets."""
    species = ['Fe', 'Fe']
    magmom = [2.0, (2.0, 0.0, 0.0)]
    new_species, mapping = create_additional_species(species, magmom)
    assert new_species == ['Fe1', 'Fe2']
    assert mapping == {'Fe1': 2.0, 'Fe2': (2.0, 0.0, 0.0)}


def test_create_additional_species_mismatched_length():
    """Mismatch between species and magmom should raise."""
    import pytest

    with pytest.raises(ValueError):
        create_additional_species(['Fe', 'O'], [1.0])


def test_convert_to_plain_list_roundtrip_scalar():
    species = ['Fe', 'O', 'Fe', 'O']
    magmom = [2.0, 0.0, -2.0, 0.0]
    new_species, mapping = create_additional_species(species, magmom)
    symbols, magmoms = convert_to_plain_list(new_species, mapping)
    assert symbols == ['Fe', 'O', 'Fe', 'O']
    assert magmoms == [2.0, 0.0, -2.0, 0.0]


def test_convert_to_plain_list_roundtrip_vector():
    species = ['Fe', 'O', 'Fe', 'O']
    magmom = [(2.0, 0.0, 0.0), (0.0, 0.0, 0.0), (-2.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
    new_species, mapping = create_additional_species(species, magmom)
    symbols, magmoms = convert_to_plain_list(new_species, mapping)
    assert symbols == ['Fe', 'O', 'Fe', 'O']
    assert magmoms == [(2.0, 0.0, 0.0), (0.0, 0.0, 0.0), (-2.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
