"""
Default input sets for VASP
"""

from copy import deepcopy
from typing import Union

from aiida.orm import Dict, StructureData

from aiida_vasp.parsers.content_parsers.potcar import MultiPotcarIo

from .base import FELEMS, InputSet


class VASPInputSet(InputSet):
    """Input set for VASP"""

    def get_input_dict(self, structure: StructureData, raw_python: bool = True) -> Union[dict, Dict]:
        """
        Compose the Dict object containing the input settings.
        """
        out_dict = super().get_input_dict(structure, raw_python=True)

        # Check if there is any magnetic elements
        spin = False
        mapping = deepcopy(self._presets['magmom_mapping'])
        # Update with overrides
        mapping.update(self.overrides.get('magmom_mapping', {}))
        default = mapping['default']
        kind_symbols = [kind.name for kind in structure.kinds]
        for symbol in mapping:
            if symbol in kind_symbols:
                spin = True
                break
        if 'magmom_mapping' in self.overrides or 'magmom' in self.overrides:
            spin = True

        # Setup magnetic moments
        magmom = []
        if spin:
            if isinstance(structure, StructureData):
                for site in structure.sites:
                    magmom.append(mapping.get(site.kind_name, default))
            else:
                for atom in structure:
                    magmom.append(mapping.get(atom.symbol, default))
        if magmom:
            out_dict['ispin'] = 2
            out_dict['magmom'] = magmom

        # Setup LDAU parameters
        ldauumap = deepcopy(self._presets['ldauu_mapping'])
        ldauumap.update(self.overrides.get('ldauu_mapping', {}))

        ldaujmap = deepcopy(self._presets['ldauj_mapping'])
        ldaujmap.update(self.overrides.get('ldauj_mapping', {}))

        ldaukeys = get_ldau_keys(structure, ldauumap, utype=2, jmapping=ldaujmap)
        out_dict.update(ldaukeys)

        # Apply overrides again over the automatically applied keys
        self.apply_overrides(out_dict)

        if not raw_python:
            out_dict = Dict(dict=out_dict)
        return out_dict

    def get_pp_mapping(self, structure: StructureData) -> dict:
        """Return the mapping from element to the POTCAR name"""
        elms = [kind.name for kind in structure.kinds]

        pmap = deepcopy(self._presets['potcar_mapping'])
        # Update the mapping from override, if any
        pmap.update(self.overrides.get('potcar_mapping', {}))

        out_dict = {key: pmap[key] for key in elms}
        return out_dict

    def get_potcar_family(self) -> str:
        return self._presets['potcar_family']

    def get_kpoints_spacing(self) -> float:
        return self._presets.get('kpoints_spacing')


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
