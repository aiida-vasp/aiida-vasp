"""
Using SUMO interface for getting the band structure
This module requires sumo to be installed.
"""
from sumo.symmetry.kpoints import get_path_data

from aiida import orm
from aiida.engine import calcfunction


@calcfunction
def kpath_from_sumo(structure: orm.StructureData, mode: orm.Str, symprec: orm.Float, line_density):
    """
    Obtain kpoint path from sumo

    Supports multiple modes: bradcrack, pymatgen, latimer-munro, seekpath
    """

    struct = structure.get_pymatgen()
    line_density = line_density.value

    path, kpoints_raw, labels = get_path_data(
        struct,
        mode.value,
        symprec.value,
        line_density=line_density,
    )
    # Primitive structure
    prim = orm.StructureData(pymatgen=path.prim)

    # kpoints
    kpoints = orm.KpointsData()
    kpoints.set_kpoints(kpoints_raw)

    actual_labels = []
    for idx, label in enumerate(labels):
        if label != '':
            # Standarise GAMMA handling
            if 'GAMMA' in label:
                label = 'GAMMA'
            actual_labels.append([idx, label])
    # Set label locations
    kpoints.labels = actual_labels

    return {'primitive_structure': prim, 'explicit_kpoints': kpoints}
