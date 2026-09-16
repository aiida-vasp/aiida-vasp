"""Unit tests for the MAGMOM helper functions used by the v2 VaspWorkChain."""

import pytest

from aiida_vasp.workchains.v2.vasp import (
    _is_vector_magmom,
    _magmom_entry_to_components,
    _magmom_to_incar,
    _normalise_per_atom_magmom,
    _set_spin_parameters,
)


def test_is_vector_magmom():
    """A magmom entry is a 3-vector iff it is a list/tuple of length 3."""
    assert _is_vector_magmom([2.0, 0.0, 0.0]) is True
    assert _is_vector_magmom((2.0, 0.0, 0.0)) is True
    assert _is_vector_magmom([2.0, 0.0]) is False
    assert _is_vector_magmom([2.0, 0.0, 0.0, 1.0]) is False
    assert _is_vector_magmom(2.0) is False
    assert _is_vector_magmom('2.0') is False


def test_magmom_entry_to_components():
    """Scalar entries become one component; 3-vectors become three floats."""
    assert _magmom_entry_to_components(2.0) == [2.0]
    assert _magmom_entry_to_components(2) == [2.0]
    assert _magmom_entry_to_components([2.0, 0.0, 0.0]) == [2.0, 0.0, 0.0]
    assert _magmom_entry_to_components((2.0, 0.0, 0.0)) == [2.0, 0.0, 0.0]
    with pytest.raises(ValueError, match='exactly 3 components'):
        _magmom_entry_to_components([2.0, 0.0])


def test_normalise_per_atom_magmom_nested():
    """Nested list-of-lists is returned as per-site 3-vectors."""
    result = _normalise_per_atom_magmom([[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]])
    assert result == [[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]]


def test_normalise_per_atom_magmom_flat_with_n_sites():
    """Flat 3-vector representation is regrouped when the site count is known."""
    result = _normalise_per_atom_magmom([2.0, 0.0, 0.0, -2.0, 0.0, 0.0], n_sites=2)
    assert result == [[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]]


def test_normalise_per_atom_magmom_scalar_and_empty():
    """Scalar per-site entries pass through unchanged; empty input stays empty."""
    assert _normalise_per_atom_magmom([2.0, -2.0]) == [2.0, -2.0]
    assert _normalise_per_atom_magmom([]) == []
    # Without a known site count a flat list is returned verbatim
    assert _normalise_per_atom_magmom([2.0, 0.0, 0.0]) == [2.0, 0.0, 0.0]


def test_magmom_to_incar():
    """Magmom lists are serialised to a single space-separated INCAR string."""
    assert _magmom_to_incar([2.0, -2.0]) == '2.0 -2.0'
    assert _magmom_to_incar([(2.0, 0.0, 0.0), (-2.0, 0.0, 0.0)]) == '2.0 0.0 0.0 -2.0 0.0 0.0'
    assert _magmom_to_incar([[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]]) == '2.0 0.0 0.0 -2.0 0.0 0.0'
    with pytest.raises(ValueError, match='all scalars or all 3-component vectors'):
        _magmom_to_incar([2.0, (1.0, 0.0, 0.0)])
    # Already serialised strings pass through untouched
    assert _magmom_to_incar('2.0 -2.0') == '2.0 -2.0'


def test_set_spin_parameters():
    """Scalar and vector moments enable their corresponding spin modes."""
    scalar_parameters = {}
    _set_spin_parameters(scalar_parameters, [2.0, -2.0])
    assert scalar_parameters == {'ispin': 2}

    explicit_ispin = {'ispin': 1}
    _set_spin_parameters(explicit_ispin, [2.0, -2.0])
    assert explicit_ispin == {'ispin': 1}

    vector_parameters = {}
    _set_spin_parameters(vector_parameters, [(2.0, 0.0, 0.0), (-2.0, 0.0, 0.0)])
    assert vector_parameters == {'lnoncollinear': True}

    soc_parameters = {'lsorbit': True}
    _set_spin_parameters(soc_parameters, [(2.0, 0.0, 0.0)])
    assert soc_parameters == {'lsorbit': True}
