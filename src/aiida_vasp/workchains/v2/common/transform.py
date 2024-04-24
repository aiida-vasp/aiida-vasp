"""
Collection of process functions for AiiDA, used for structure transformation
"""
import re
from typing import List, Tuple

from ase import Atoms
import numpy as np

from aiida import orm
from aiida.engine import calcfunction
from aiida.orm import ArrayData, CalcFunctionNode, Node, QueryBuilder, StructureData


@calcfunction
def magnetic_structure_decorate(structure, magmom):
    """
    Create Quantum Espresso style decorated structure with
    given magnetic moments.
    """

    magmom = magmom.get_list()
    assert (
        len(magmom) == len(structure.sites)
    ), f"Mismatch between the magmom ({len(magmom)}) and the nubmer of sites ({len(structure.sites)})."
    old_species = [structure.get_kind(site.kind_name).symbol for site in structure.sites]
    new_species, magmom_mapping = create_additional_species(old_species, magmom)
    new_structure = StructureData()
    new_structure.set_cell(structure.cell)
    new_structure.set_pbc(structure.pbc)
    for site, name in zip(structure.sites, new_species):
        this_symbol = structure.get_kind(site.kind_name).symbol
        new_structure.append_atom(position=site.position, symbols=this_symbol, name=name)

    # Keep the label
    new_structure.label = structure.label
    return {'structure': new_structure, 'mapping': orm.Dict(dict=magmom_mapping)}


@calcfunction
def magnetic_structure_dedecorate(structure, mapping):
    """
    Remove decorations of a structure with multiple names for the same specie
    given that the decoration was previously created to give different species
    name for different initialisation of magnetic moments.
    """

    mapping = mapping.get_dict()
    # Get a list of decroated names
    old_species = [structure.get_kind(site.kind_name).name for site in structure.sites]
    new_species, magmom = convert_to_plain_list(old_species, mapping)

    new_structure = StructureData()
    new_structure.set_cell(structure.cell)
    new_structure.set_pbc(structure.pbc)

    for site, name in zip(structure.sites, new_species):
        this_symbol = structure.get_kind(site.kind_name).symbol
        new_structure.append_atom(position=site.position, symbols=this_symbol, name=name)
    new_structure.label = structure.label
    return {'structure': new_structure, 'magmom': orm.List(list=magmom)}


@calcfunction
def rattle(structure, amp):
    """
    Rattle the structure by a certain amplitude
    """
    native_keys = ['cell', 'pbc1', 'pbc2', 'pbc3', 'kinds', 'sites', 'mp_id']
    # Keep the foreign keys as it is
    foreign_attrs = {key: value for key, value in structure.attributes.items() if key not in native_keys}
    atoms = structure.get_ase()
    atoms.rattle(amp.value)
    # Clean any tags etc
    atoms.set_tags(None)
    atoms.set_masses(None)
    # Convert it back
    out = StructureData(ase=atoms)
    out.base.attributes.set_many(foreign_attrs)
    out.label = structure.label + ' RATTLED'
    return out


@calcfunction
def get_primitive(structure):
    """Create primitive structure use pymatgen interface"""
    from aiida.orm import StructureData

    pstruct = structure.get_pymatgen()
    ps = pstruct.get_primitive_structure()
    out = StructureData(pymatgen=ps)
    out.label = structure.label + ' PRIMITIVE'
    return out


@calcfunction
def get_standard_primitive(structure, **kwargs):
    """Create the standard primitive structure via seekpath"""
    from aiida.tools.data.array.kpoints import get_kpoints_path

    parameters = kwargs.get('parameters', {'symprec': 1e-5})
    if isinstance(parameters, orm.Dict):
        parameters = parameters.get_dict()

    out = get_kpoints_path(structure, **parameters)['primitive_structure']
    out.label = structure.label + ' PRIMITIVE'
    return out


@calcfunction
def spglib_refine_cell(structure, symprec):
    """Create the standard primitive structure via seekpath"""
    from spglib import refine_cell

    from aiida.tools.data.structure import spglib_tuple_to_structure, structure_to_spglib_tuple

    structure_tuple, kind_info, kinds = structure_to_spglib_tuple(structure)

    lattice, positions, types = refine_cell(structure_tuple, symprec.value)

    refined = spglib_tuple_to_structure((lattice, positions, types), kind_info, kinds)

    return refined


@calcfunction
def get_standard_conventional(structure):
    """Create the standard primitive structure via seekpath"""
    from aiida.tools.data.array.kpoints import get_kpoints_path

    out = get_kpoints_path(structure)['conv_structure']
    out.label = structure.label + ' PRIMITIVE'
    return out


@calcfunction
def get_refined_structure(structure, symprec, angle_tolerance):
    """Create refined structure use pymatgen's interface"""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    from aiida.orm import StructureData

    pstruct = structure.get_pymatgen()
    ana = SpacegroupAnalyzer(pstruct, symprec=symprec.value, angle_tolerance=angle_tolerance.value)
    ps = ana.get_refined_structure()
    out = StructureData(pymatgen=ps)
    out.label = structure.label + ' REFINED'
    return out


@calcfunction
def get_conventional_standard_structure(structure, symprec, angle_tolerance):
    """Create conventional standard structure use pymatgen's interface"""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    from aiida.orm import StructureData

    pstruct = structure.get_pymatgen()
    ana = SpacegroupAnalyzer(pstruct, symprec=symprec.value, angle_tolerance=angle_tolerance.value)
    ps = ana.get_conventional_standard_structure()
    out = StructureData(pymatgen=ps)
    out.label = structure.label + ' CONVENTIONAL STANDARD'
    return out


@calcfunction
def make_supercell(structure, supercell, **kwargs):
    """Make supercell structure, keep the tags in order"""
    from ase.build.supercells import make_supercell as ase_supercell

    if 'tags' in kwargs:
        tags = kwargs['tags']
    else:
        tags = None

    atoms = structure.get_ase()
    atoms.set_tags(tags)

    slist = supercell.get_list()
    if isinstance(slist[0], int):
        satoms = atoms.repeat(slist)
    else:
        satoms = ase_supercell(atoms, np.array(slist))
    if 'no_sort' not in kwargs:
        satoms = sort(satoms)

    if tags:
        stags = satoms.get_tags().tolist()
    satoms.set_tags(None)

    out = StructureData(ase=satoms)
    out.label = structure.label + ' SUPER {} {} {}'.format(*slist)

    if tags:
        return {'structure': out, 'tags': orm.List(list=stags)}
    else:
        return {'structure': out}


@calcfunction
def niggli_reduce(structure):
    """Peroform niggli reduction using ase as the backend - this will rotate the structure into the standard setting"""
    from ase.build import niggli_reduce as niggli_reduce_

    atoms = structure.get_ase()
    niggli_reduce_(atoms)
    new_structure = StructureData(ase=atoms)
    new_structure.label = structure.label + ' NIGGLI'
    return new_structure


@calcfunction
def niggli_reduce_spglib(structure):
    """Peroform niggli reduction using spglib as backend - this does not rotate the structure"""
    from spglib import niggli_reduce as niggli_reduce_spg

    atoms = structure.get_ase()
    reduced_cell = niggli_reduce_spg(atoms.cell)
    atoms.set_cell(reduced_cell)
    atoms.wrap()
    new_structure = StructureData(ase=atoms)
    new_structure.label = structure.label + ' NIGGLI'
    return new_structure


@calcfunction
def neb_interpolate(init_structure, final_strucrture, nimages):
    """
    Interplate NEB frames using the starting and the final structures

    Get around the PBC warpping problem by calculating the MIC displacements
    from the initial to the final structure
    """
    from ase.neb import NEB

    ainit = init_structure.get_ase()
    afinal = final_strucrture.get_ase()
    disps = []

    # Find distances
    acombined = ainit.copy()
    acombined.extend(afinal)
    # Get piece-wise MIC distances
    for i in range(len(ainit)):
        dist = acombined.get_distance(i, i + len(ainit), vector=True, mic=True)
        disps.append(dist.tolist())
    disps = np.asarray(disps)
    ainit.wrap(eps=1e-1)
    afinal = ainit.copy()

    # Displace the atoms according to MIC distances
    afinal.positions += disps
    neb = NEB([ainit.copy() for i in range(int(nimages) + 1)] + [afinal.copy()])
    neb.interpolate()
    out_init = StructureData(ase=neb.images[0])
    out_init.label = init_structure.label + ' INIT'
    out_final = StructureData(ase=neb.images[-1])
    out_final.label = init_structure.label + ' FINAL'

    outputs = {'image_init': out_init}
    for i, out in enumerate(neb.images[1:-1]):
        outputs[f"image_{i+1:02d}"] = StructureData(ase=out)
        outputs[f"image_{i+1:02d}"].label = init_structure.label + f" FRAME {i+1:02d}"
    outputs['image_final'] = out_final
    return outputs


@calcfunction
def fix_atom_order(reference, to_fix):
    """
    Fix atom order by finding NN distances bet ween two frames. This resolves
    the issue where two closely matching structures having diffferent atomic orders.
    Note that the two frames must be close enough for this to work
    """

    aref = reference.get_ase()
    afix = to_fix.get_ase()

    # Index of the reference atom in the second structure
    new_indices = np.zeros(len(aref), dtype=int)

    # Find distances
    acombined = aref.copy()
    acombined.extend(afix)
    # Get piece-wise MIC distances
    for i in range(len(aref)):
        dists = []
        for j in range(len(aref)):
            dist = acombined.get_distance(i, j + len(aref), mic=True)
            dists.append(dist)
        min_idx = np.argmin(dists)
        min_dist = min(dists)
        if min_dist > 0.5:
            print(f"Large displacement found - moving atom {j} to {i} - please check if this is correct!")
        new_indices[i] = min_idx

    afixed = afix[new_indices]
    fixed_structure = StructureData(ase=afixed)
    fixed_structure.label = to_fix.label + ' UPDATED ORDER'
    return fixed_structure


def match_atomic_order_(atoms: Atoms, atoms_ref: Atoms) -> Tuple[Atoms, List[int]]:
    """
    Reorder the atoms to that of the reference.

    Only works for identical or nearly identical structures that are ordered differently.
    Returns a new `Atoms` object with order similar to that of `atoms_ref` as well as the sorting indices.
    """

    # Find distances
    acombined = atoms_ref.copy()
    acombined.extend(atoms)
    new_index = []
    # Get piece-wise MIC distances
    jidx = list(range(len(atoms), len(atoms) * 2))
    for i in range(len(atoms)):
        dists = acombined.get_distances(i, jidx, mic=True)
        # Find the index of the atom with the smallest distance
        min_idx = np.where(dists == dists.min())[0][0]
        new_index.append(min_idx)
    assert len(set(new_index)) == len(atoms), 'The detected mapping is not unique!'
    return atoms[new_index], new_index


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
    for symbol, magmom in zip(species, magmom):
        current_symbol = symbol
        # Mappings for this original symbol
        mapping = current_species_mapping[symbol]
        # First check if this magmom has been treated
        not_seen = True
        for sym_, mag_ in mapping.items():
            if mag_ == magmom:
                current_symbol = sym_
                not_seen = False
        # This symbol has not been seen yet
        if not_seen:
            if current_symbol in mapping:
                # The other species having the same symbol has been assigned
                counter = len(mapping) + 1
                current_symbol = f"{symbol}{counter}"
            mapping[current_symbol] = magmom
        new_species.append(current_symbol)

    # Rename symbols that has more than one species, so A becomes A1
    for symbol, mapping in current_species_mapping.items():
        if len(mapping) > 1:
            mapping[f"{symbol}1"] = mapping[symbol]
            mapping.pop(symbol)
            # Refresh the new_species list
            new_species = [f"{sym}1" if sym == symbol else sym for sym in new_species]

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
        match = re.match(r'(\w+)\d+', symbol)
        if match:
            symbol = match.group(1)
        symbols.append(symbol)
    return symbols, magmoms
