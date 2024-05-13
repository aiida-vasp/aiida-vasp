"""
Module for converting between different representation for
initial site magnetization.
"""

import re


def create_additional_species(species: list, magmom: list):
    """
    Create additional species depending on magnetic moments.
    For example, create Fe1 and Fe2 if there are Fe with different
    magnetisations.

    Returns:
        a tuples of (newspecies, magmom_mapping)
    """

    unique_species = set(species)
    new_species = []
    current_species_mapping = {sym: {} for sym in unique_species}
    for symbol, this_mag in zip(species, magmom):
        current_symbol = symbol
        # Mappings for this original symbol
        mapping = current_species_mapping[symbol]
        # First check if this magmom has been treated
        not_seen = True
        for sym_, mag_ in mapping.items():
            if mag_ == this_mag:
                current_symbol = sym_
                not_seen = False
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

    all_mapping = {}
    for value in current_species_mapping.values():
        all_mapping.update(value)

    return new_species, all_mapping


def convert_to_plain_list(species: list, magmom_mapping: dict):
    """
    Covert from a decorated species list to a plain list of symbols
    and magnetic moments.

    Returns:
        A tuple of (symbols, magmoms)
    """
    magmoms = []
    symbols = []
    for symbol in species:
        magmoms.append(magmom_mapping[symbol])
        # Drop the number suffix in the symbol
        match = re.match(r'(\w+)\d+', symbol)
        if match:
            symbols.append(match.group(1))
        else:
            symbols.append(symbol)
    return symbols, magmoms
