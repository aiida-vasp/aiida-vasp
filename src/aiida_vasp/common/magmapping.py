"""
Module for converting between different representation for
initial site magnetization.
"""

from __future__ import annotations

import re
from typing import List, Sequence, Tuple, Union

# A magmom can be either a scalar (collinear) or a 3-component vector
# (non-collinear / spin-orbit). We accept anything that can be turned into a
# tuple of floats by ``_normalise_magmom``.
MagmomScalar = float
MagmomVector = Tuple[float, float, float]
Magmom = Union[MagmomScalar, MagmomVector, Sequence[float]]


def _normalise_magmom(value: Magmom) -> Union[float, Tuple[float, float, float]]:
    """Convert a magmom into a canonical form.

    A scalar is kept as a float. A sequence is converted into a 3-tuple (any
    trailing components are ignored, missing components default to 0.0).
    """
    if isinstance(value, (int, float)):
        return float(value)
    seq = list(value)
    if len(seq) == 0:
        raise ValueError('Magmom must have at least one component.')
    if len(seq) == 1:
        return float(seq[0])
    # Pad/truncate to 3 components
    seq = list(seq[:3]) + [0.0] * max(0, 3 - len(seq))
    return (float(seq[0]), float(seq[1]), float(seq[2]))


def _magmom_equal(a: Union[float, Tuple[float, float, float]], b: Union[float, Tuple[float, float, float]]) -> bool:
    """Equality check for normalised magmoms."""
    if isinstance(a, float) and isinstance(b, float):
        return a == b
    if isinstance(a, tuple) and isinstance(b, tuple):
        return a == b
    # Mixed scalar vs vector: never equal (different species buckets)
    return False


def create_additional_species(
    species: List[str], magmom: List[Magmom]
) -> Tuple[List[str], dict[str, Union[float, Tuple[float, float, float]]]]:
    """
    Create additional species depending on magnetic moments.

    For example, create ``Fe1`` and ``Fe2`` if there are Fe atoms with
    different magnetic moments. Works for both scalar magmoms (collinear
    case) and 3-component magmom vectors (non-collinear / spin-orbit case
    with ``LNONCOLLINEAR = .TRUE.``).

    :param species: list of element symbols, one per site.
    :param magmom: list of magnetic moments, one per site. Each entry may be a
        scalar or a 3-component sequence (x, y, z).
    :return: a tuple of (new_species, magmom_mapping). The mapping contains one
        entry per generated species, mapping the decorated name to the
        normalised magmom (float or 3-tuple).
    """

    if len(species) != len(magmom):
        raise ValueError(
            f'species and magmom must have the same length, got {len(species)} and {len(magmom)}.'
        )

    unique_species = set(species)
    new_species: List[str] = []
    current_species_mapping: dict[str, dict[str, Union[float, Tuple[float, float, float]]]] = {
        sym: {} for sym in unique_species
    }
    for symbol, raw_mag in zip(species, magmom):
        this_mag = _normalise_magmom(raw_mag)
        current_symbol = symbol
        # Mappings for this original symbol
        mapping = current_species_mapping[symbol]
        # First check if this magmom has been treated
        not_seen = True
        for sym_, mag_ in mapping.items():
            if _magmom_equal(mag_, this_mag):
                current_symbol = sym_
                not_seen = False
                break
        # This symbol has not been seen yet
        if not_seen:
            if current_symbol in mapping:
                # The other species having the same symbol has been assigned
                counter = len(mapping) + 1
                current_symbol = f'{symbol}{counter}'
            mapping[current_symbol] = this_mag
        new_species.append(current_symbol)

    # Rename symbols that has more than one species, so A becomes A1
    for symbol, mapping in current_species_mapping.items():
        if len(mapping) > 1:
            mapping[f'{symbol}1'] = mapping[symbol]
            mapping.pop(symbol)
            # Refresh the new_species list
            new_species = [f'{sym}1' if sym == symbol else sym for sym in new_species]

    all_mapping: dict[str, Union[float, Tuple[float, float, float]]] = {}
    for value in current_species_mapping.values():
        all_mapping.update(value)

    return new_species, all_mapping


def convert_to_plain_list(
    species: List[str], magmom_mapping: dict[str, Magmom]
) -> Tuple[List[str], List[Union[float, Tuple[float, float, float]]]]:
    """
    Convert from a decorated species list to a plain list of symbols and
    magnetic moments.

    :return: a tuple of (symbols, magmoms).
    """
    magmoms: List[Union[float, Tuple[float, float, float]]] = []
    symbols: List[str] = []
    for symbol in species:
        if symbol not in magmom_mapping:
            raise ValueError(f'No magmom mapping entry for species {symbol!r}.')
        magmoms.append(_normalise_magmom(magmom_mapping[symbol]))
        # Drop the number suffix in the symbol
        match = re.match(r'(\w+)\d+', symbol)
        if match:
            symbols.append(match.group(1))
        else:
            symbols.append(symbol)
    return symbols, magmoms
