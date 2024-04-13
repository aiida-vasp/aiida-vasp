"""
Utility functions for running NEB calculations
"""
from ase.neb import NEB
import numpy as np

from aiida.engine import calcfunction
from aiida.orm import StructureData


@calcfunction
def neb_interpolate(init_structure, final_strucrture, nimages):
    """
    Interplate NEB frames using the starting and the final structures

    Get around the PBC warpping problem by calculating the MIC displacements
    from the initial to the final structure
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
        outputs[f'image_{i+1:02d}'] = StructureData(ase=out)
        outputs[f'image_{i+1:02d}'].label = init_structure.label + f' FRAME {i+1:02d}'
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
            print(f'Large displacement found - moving atom {j} to {i} - please check if this is correct!')
        new_indices[i] = min_idx

    afixed = afix[new_indices]
    fixed_structure = StructureData(ase=afixed)
    fixed_structure.label = to_fix.label + ' UPDATED ORDER'
    return fixed_structure
