"""
Collection of process functions for AiiDA, used for structure transformation
"""

from typing import Any, List, Tuple, Union

import numpy as np
from aiida import orm
from aiida.engine import calcfunction
from aiida.orm import StructureData
from aiida.tools.data.array.kpoints import get_kpoints_path
from aiida.tools.data.structure import spglib_tuple_to_structure, structure_to_spglib_tuple
from ase import Atoms
from ase.build import niggli_reduce as niggli_reduce_
from ase.build import sort
from ase.build.supercells import make_supercell as ase_supercell

try:
    from ase.mep.neb import NEB
except ModuleNotFoundError:
    from ase.neb import NEB

from spglib import niggli_reduce as niggli_reduce_spg
from spglib import refine_cell

from aiida_vasp.common.magmapping import convert_to_plain_list, create_additional_species


def _coerce_magmom_list(raw: list[Any]) -> list[Any]:
    """Coerce an ``orm.List`` payload into the canonical magmom representation.

    AiiDA flattens lists-of-lists into flat lists when stored in an
    :class:`orm.List`. This helper accepts both forms (``[[1.0, 0.0, 0.0], 2.0]``
    or ``[1.0, 0.0, 0.0, 2.0]``) and re-groups the values per site.
    """
    if not raw:
        return list(raw)

    # Detect if a nested structure was preserved (each entry is itself a list).
    if all(isinstance(entry, (list, tuple)) for entry in raw):
        coerced: list[Any] = []
        for entry in raw:
            coerced.append(_normalise_entry(entry))
        return coerced

    # Otherwise, try to detect flat 3-vectors. A flat list will have a length
    # that is a multiple of the number of sites. We treat the case where the
    # length is exactly ``3 * n_sites`` as a flat 3-vector representation.
    # Without additional context (the structure), we cannot decide here, so we
    # return the flat list as-is; downstream code should provide a nested
    # representation for 3D magmoms.
    return list(raw)


def _normalise_entry(entry: Any) -> Union[float, Tuple[float, float, float]]:
    """Normalise a single magmom entry."""
    if isinstance(entry, (int, float)):
        return float(entry)
    seq = list(entry)
    if len(seq) == 0:
        raise ValueError('Magmom must have at least one component.')
    if len(seq) == 1:
        return float(seq[0])
    seq = list(seq[:3]) + [0.0] * max(0, 3 - len(seq))
    return (float(seq[0]), float(seq[1]), float(seq[2]))


@calcfunction
def magnetic_structure_decorate(structure: orm.StructureData, magmom: orm.List) -> dict[str, Any]:
    """
    Create Quantum Espresso style decorated structure with
    given magnetic moments.

    Each entry of ``magmom`` may be a scalar (collinear) or a 3-component
    sequence for ``LNONCOLLINEAR = .TRUE.``.
    """

    raw_magmom = magmom.get_list()
    magmom_list = _coerce_magmom_list(raw_magmom)
    assert len(magmom_list) == len(structure.sites), (
        f'Mismatch between the magmom ({len(magmom_list)}) and the nubmer of sites ({len(structure.sites)}).'
    )
    old_species = [structure.get_kind(site.kind_name).symbol for site in structure.sites]
    new_species, magmom_mapping = create_additional_species(old_species, magmom_list)
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
def magnetic_structure_dedecorate(structure: orm.StructureData, mapping: orm.Dict) -> dict[str, Any]:
    """
    Remove decorations of a structure with multiple names for the same specie
    given that the decoration was previously created to give different species
    name for different initialisation of magnetic moments.
    """

    raw_mapping = mapping.get_dict()
    # The mapping values are stored as plain python objects (floats or lists
    # of 3 floats). Re-normalise them so downstream code gets a consistent
    # representation.
    magmom_mapping = {key: _normalise_entry(value) for key, value in raw_mapping.items()}
    # Get a list of decroated names
    old_species = [structure.get_kind(site.kind_name).name for site in structure.sites]
    new_species, magmom = convert_to_plain_list(old_species, magmom_mapping)

    new_structure = StructureData()
    new_structure.set_cell(structure.cell)
    new_structure.set_pbc(structure.pbc)

    for site, name in zip(structure.sites, new_species):
        this_symbol = structure.get_kind(site.kind_name).symbol
        new_structure.append_atom(position=site.position, symbols=this_symbol, name=name)
    new_structure.label = structure.label
    # Store as nested list so 3D magmoms survive AiiDA's flattening in List nodes.
    return {'structure': new_structure, 'magmom': orm.List(list=list(magmom))}


@calcfunction
def rattle(structure: orm.StructureData, amp: orm.Float) -> orm.StructureData:
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
def random_rattle(structure: orm.StructureData, amp: orm.Float) -> orm.StructureData:
    """
    Rattle the structure by a certain amplitude

    This function is similar to `rattle`, but it uses a random seed so it generates
    different structures when called with the same structure.
    """
    native_keys = ['cell', 'pbc1', 'pbc2', 'pbc3', 'kinds', 'sites', 'mp_id']
    # Keep the foreign keys as it is
    foreign_attrs = {key: value for key, value in structure.attributes.items() if key not in native_keys}
    atoms = structure.get_ase()
    atoms.rattle(amp.value, np.random.randint(0, 100000000))
    # Clean any tags etc
    atoms.set_tags(None)
    atoms.set_masses(None)
    # Convert it back
    out = StructureData(ase=atoms)
    out.base.attributes.set_many(foreign_attrs)
    out.label = structure.label + ' RATTLED'
    return out


@calcfunction
def get_primitive(structure: orm.StructureData) -> orm.StructureData:
    """Create primitive structure use pymatgen interface"""
    pstruct = structure.get_pymatgen()
    ps = pstruct.get_primitive_structure()
    out = StructureData(pymatgen=ps)
    out.label = structure.label + ' PRIMITIVE'
    return out


@calcfunction
def get_standard_primitive(structure: orm.StructureData, **kwargs: Any) -> orm.StructureData:
    """Create the standard primitive structure via seekpath"""

    parameters = kwargs.get('parameters', {'symprec': 1e-5})
    if isinstance(parameters, orm.Dict):
        parameters = parameters.get_dict()

    out = get_kpoints_path(structure, **parameters)['primitive_structure']
    out.label = structure.label + ' PRIMITIVE'
    return out


@calcfunction
def spglib_refine_cell(structure: orm.StructureData, symprec: Any) -> orm.StructureData:
    """Create the standard primitive structure via seekpath"""

    structure_tuple, kind_info, kinds = structure_to_spglib_tuple(structure)

    lattice, positions, types = refine_cell(structure_tuple, symprec.value)

    refined = spglib_tuple_to_structure((lattice, positions, types), kind_info, kinds)

    return refined


@calcfunction
def get_standard_conventional(structure: orm.StructureData) -> orm.StructureData:
    """Create the standard primitive structure via seekpath"""
    out = get_kpoints_path(structure)['conv_structure']
    out.label = structure.label + ' PRIMITIVE'
    return out


@calcfunction
def get_refined_structure(structure: orm.StructureData, symprec: orm.Float, angle_tolerance: Any) -> orm.StructureData:
    """Create refined structure use pymatgen's interface"""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer  # noqa: PLC0415

    pstruct = structure.get_pymatgen()
    ana = SpacegroupAnalyzer(pstruct, symprec=symprec.value, angle_tolerance=angle_tolerance.value)
    ps = ana.get_refined_structure()
    out = StructureData(pymatgen=ps)
    out.label = structure.label + ' REFINED'
    return out


@calcfunction
def get_conventional_standard_structure(structure: orm.StructureData, symprec: orm.Float, angle_tolerance: Any) -> Any:
    """Create conventional standard structure use pymatgen's interface"""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer  # noqa: PLC0415

    pstruct = structure.get_pymatgen()
    ana = SpacegroupAnalyzer(pstruct, symprec=symprec.value, angle_tolerance=angle_tolerance.value)
    ps = ana.get_conventional_standard_structure()
    out = StructureData(pymatgen=ps)
    out.label = structure.label + ' CONVENTIONAL STANDARD'
    return out


@calcfunction
def make_supercell(structure: Any, supercell: Any, **kwargs: Any) -> dict[str, Any]:
    """Make supercell structure, keep the tags in order"""

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
def niggli_reduce(structure: orm.StructureData) -> orm.StructureData:
    """Peroform niggli reduction using ase as the backend - this will rotate the structure into the standard setting"""
    atoms = structure.get_ase()
    niggli_reduce_(atoms)
    new_structure = StructureData(ase=atoms)
    new_structure.label = structure.label + ' NIGGLI'
    return new_structure


@calcfunction
def niggli_reduce_spglib(structure: orm.StructureData) -> orm.StructureData:
    """Peroform niggli reduction using spglib as backend - this does not rotate the structure"""
    atoms = structure.get_ase()
    reduced_cell = niggli_reduce_spg(atoms.cell)
    atoms.set_cell(reduced_cell)
    atoms.wrap()
    new_structure = StructureData(ase=atoms)
    new_structure.label = structure.label + ' NIGGLI'
    return new_structure


@calcfunction
def neb_interpolate(
    init_structure: orm.StructureData, final_strucrture: orm.StructureData, nimages: orm.Int
) -> dict[str, orm.StructureData]:
    """
    Interpolate NEB frames using the starting and the final structures

    Get around the PBC wrapping problem by calculating the MIC displacements
    from the initial to the final structure

    The initial structure is not changed, while the final structure is
    modified to be consistent with the initial structure in terms of absolute
    displacements, i.e. the final structure is *unwrapped*.
    """
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
        outputs[f'image_{i + 1:02d}'] = StructureData(ase=out)
        outputs[f'image_{i + 1:02d}'].label = init_structure.label + f' FRAME {i + 1:02d}'
    outputs['image_final'] = out_final
    return outputs


@calcfunction
def fix_atom_order(reference: orm.StructureData, to_fix: orm.StructureData) -> orm.StructureData:
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
            print(f'Large displacement found - moving atom {j} to {i} - please check if this is correct!')
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
