"""
Functions to setting up LDA+U calculations
"""

from aiida.orm import StructureData

from aiida_vasp.parsers.content_parsers.potcar import MultiPotcarIo

FELEMS = [
    'La',
    'Ce',
    'Pr',
    'Nd',
    'Pm',
    'Sm',
    'Eu',
    'Gd',
    'Tb',
    'Dy',
    'Ho',
    'Er',
    'Tm',
    'Yb',
    'Lu',
    'Ac',
    'Th',
    'Pa',
    'U',
    'Np',
    'Pu',
    'Am',
    'Cm',
    'Bk',
    'Cf',
    'Es',
    'Fm',
    'Md',
    'No',
    'Lr',
]


def get_ldau_keys(
    structure: StructureData,
    mapping: dict[str, list[int | float | str]],
    utype: int = 2,
    jmapping: dict[str, list[int | float | str]] | None = None,
    felec: bool = False,
) -> dict[str, str | float]:
    """
    Setup LDAU mapping. In VASP, the U for each species has to be
    defined in the order that they appear in POSCAR. This is a helper
    to make sure the values of U are associated to each specie

    Arguments:
        structure: the structure, either StructureData or ase.Atoms is fine
        mapping: a dictionary in the format of  {"Mn": [d, 4]...} for U
        utype: the type of LDA+U, default to 2, which is the one with only one parameter
        jmapping: a dictionary in the format of  {"Mn": [d, 4]...} but for J
        felec: Whether we are dealing with f electrons, will increase lmaxmix if we are.


    Returns:
        dict_update: a dictionary to be used to update the raw input parameters for VASP
    """
    if isinstance(structure, StructureData):
        species = MultiPotcarIo.potentials_order(structure)
    else:
        # For ASE atoms, we keep the order of species occurrence no sorting is done
        species = []
        for symbol in structure.get_chemical_symbols():
            if symbol not in species:
                species.append(symbol)

    lsymbols = {'d': 2, 'f': 3, 'p': 1}
    if jmapping is None:
        jmapping = {}
    # Setup U array
    ldauu = []
    ldauj = []
    ldaul = []
    count = 0
    for specie in species:
        if specie in mapping:
            uvalue = mapping[specie][1]
            j = jmapping.get(specie, 0.0)
            ldaul.append(lsymbols[mapping[specie][0]])
            ldauu.append(mapping[specie][1])

            j = jmapping.get(specie, 0.0)
            ldauj.append(j)

            if specie in FELEMS:
                felec = True
            # Count the number of valid mappings
            if uvalue != 0.0 or j != 0.0:
                count += 1

        else:
            ldauu.append(0.0)
            ldauj.append(0.0)
            ldaul.append(-1)

    if count > 0:
        # Only enable U is there is any non-zero value
        output = {
            'ldauu': ldauu,
            'ldauj': ldauj,
            'ldautype': utype,
            'lmaxmix': 6 if felec else 4,
            'ldaul': ldaul,
            'ldau': True,
        }
    else:
        output = {}
    return output
