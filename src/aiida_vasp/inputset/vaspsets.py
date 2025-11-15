"""
Default input sets for VASP
"""

from copy import deepcopy
from typing import Union

from aiida.orm import Dict, StructureData

from aiida_vasp.utils.ldau import get_ldau_keys

from .base import InputSet


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
