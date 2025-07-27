"""
Using SUMO interface for getting the band structure
This module requires sumo to be installed.
"""

from aiida import orm
from aiida.engine import calcfunction
from sumo.symmetry.kpoints import get_path_data


@calcfunction
def kpath_from_sumo(
    structure: orm.StructureData, mode: orm.Str, symprec: orm.Float, line_density: orm.Float
) -> dict[str, orm.StructureData | orm.KpointsData]:
    """
    Obtain kpoint path from sumo

    Supports multiple modes: bradcrack, pymatgen, latimer-munro, seekpath
    """

    struct = structure.get_pymatgen()
    line_density_value: float = line_density.value

    path, kpoints_raw, labels = get_path_data(
        struct,
        mode.value,
        symprec.value,
        line_density=line_density_value,
    )
    # Primitive structure
    prim: orm.StructureData = orm.StructureData(pymatgen=path.prim)

    # kpoints
    kpoints: orm.KpointsData = orm.KpointsData()
    kpoints.set_kpoints(kpoints_raw)

    actual_labels: list[list[int | str]] = []
    for idx, label in enumerate(labels):
        if label != '':
            # Standardise GAMMA handling
            if 'GAMMA' in label:
                actual_labels.append([idx, 'GAMMA'])
            else:
                actual_labels.append([idx, label])
    # Set label locations
    kpoints.labels = actual_labels

    return {'primitive_structure': prim, 'explicit_kpoints': kpoints}


@calcfunction
def kpath_from_sumo_v2(
    structure: orm.StructureData, band_settings: orm.Dict
) -> dict[str, orm.StructureData | orm.KpointsData]:
    """
    Obtain kpoint path from sumo

    Supports multiple modes: bradcrack, pymatgen, latimer-munro, seekpath
    """

    struct = structure.get_pymatgen()

    path, kpoints_raw, labels = get_path_data(struct, **band_settings.get_dict())
    # Primitive structure
    prim = orm.StructureData(pymatgen=path.prim)

    # kpoints
    kpoints = orm.KpointsData()
    kpoints.set_kpoints(kpoints_raw)

    actual_labels: list[list[int | str]] = []
    for idx, label in enumerate(labels):
        if label != '':
            # Standardise GAMMA handling
            if 'GAMMA' in label:
                actual_labels.append([idx, 'GAMMA'])
            else:
                actual_labels.append([idx, label])
    # Set label locations
    kpoints.labels = actual_labels

    return {'primitive_structure': prim, 'explicit_kpoints': kpoints}
