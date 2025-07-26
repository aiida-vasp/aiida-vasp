"""
Module for using pymatgen.io.vasp.sets based input sets.
"""

from typing import Dict, Optional, Union

import numpy as np
from aiida import orm

from .base import InputSet

try:
    import pymatgen.io.vasp.sets as pmg_sets
    from pymatgen.io.vasp.inputs import KpointsSupportedModes
except ImportError:
    pmg_sets = None


class PymatgenInputSet(InputSet):
    """
    Input set using pymatgen.io.vasp.sets.

    Provides basic compatibility with pymatgen sets for generating VASP input
    parameters, k-point meshes, and pseudopotential mappings.
    """

    # An none-exhaust list of supported pymatgen input sets
    KNOWN_SETS = (
        'MPRelaxSet',
        'MITRelaxSet',
        'MPScanRelaxSet',
        'MP24RelaxSet',
        'MPMetalRelaxSet',
        'MPHSERelaxSet',
        'MVLGWSet',
        'MPAbsorptionSet',
        'MatPESStaticSet',
        'MPScanStaticSet',
        'MP24StaticSet',
        'MPHSEBSSet',
        'MPNonSCFSet',
        'MPSOCSet',
        'MPNMRSet',
        'MPStaticSet',
    )

    def __init__(
        self,
        set_name: str,
        overrides: Optional[Dict] = None,
        verbose: Optional[bool] = None,
        pmg_kwargs: Optional[Dict] = None,
    ) -> None:
        """
        Instantiate a PymatgenInputSet.

        :param set_name: Name of the pymatgen input set to use
        :type set_name: str
        :param overrides: Dictionary of parameter overrides for the input set
        :type overrides: dict or None
        :param verbose: If True, print additional information during processing
        :type verbose: bool or None
        :param pmg_kwargs: Additional keyword arguments to pass to the pymatgen input set
        :type pmg_kwargs: dict or None

        :raises AssertionError: If set_name is not in KNOWN_SETS
        """
        assert set_name in self.KNOWN_SETS, f'Unsupported set name: {set_name}'
        super().__init__(set_name, overrides=overrides, verbose=verbose)
        self._pmg_kwargs = pmg_kwargs or {}

    def _load_data(self) -> None:
        """
        Load the pymatgen input set class.

        Dynamically imports and stores the pymatgen input set class based on set_name.

        :raises ImportError: If pymatgen is not installed or cannot be imported
        """
        if pmg_sets is None:
            raise ImportError('pymatgen is not installed. Please install it to use PymatgenInputSet.')
        self._pmg_class = getattr(pmg_sets, self.set_name)

    def get_input_dict(self, structure: orm.StructureData, raw_python: bool = True) -> Union[Dict, orm.Dict]:
        """
        Compute the input parameters for a VASP calculation using pymatgen.io.vasp.sets.

        Generates INCAR parameters by instantiating the pymatgen input set with the
        given structure and applying any specified overrides. Removes certain
        parameters that conflict with aiida-vasp's input validation.

        :param structure: Crystal structure for the calculation
        :type structure: orm.StructureData
        :param raw_python: If True, return a Python dict; if False, return orm.Dict
        :type raw_python: bool

        :returns: Dictionary of INCAR parameters
        :rtype: dict or orm.Dict
        """
        ps = structure.get_pymatgen()
        pmgset = self._pmg_class(ps, **self._pmg_kwargs)
        incar_dict = {key.lower(): value for key, value in pmgset.incar.items()}
        # Apply the overrides
        for key, value in self.overrides.items():
            if value is None:
                if key in incar_dict:
                    incar_dict.pop(key)
            else:
                incar_dict[key] = value

        # pop icharg which conflicts with aiida-vasp's input checks
        incar_dict.pop('icharg', None)
        incar_dict.pop('istart', None)
        incar_dict.pop('kspacing', None)

        if raw_python:
            return incar_dict
        return orm.Dict(dict=incar_dict)

    def get_pp_mapping(self, structure: orm.StructureData) -> Dict[str, str]:
        """
        Get the pseudopotential mapping used by the input set.

        Returns a dictionary mapping element symbols to their corresponding
        pseudopotential symbols as defined by the pymatgen input set.

        :param structure: Crystal structure for the calculation
        :type structure: orm.StructureData

        :returns: Dictionary mapping element names to pseudopotential symbols
        :rtype: dict
        """
        ps = structure.get_pymatgen()
        pmgset = self._pmg_class(ps, **self._pmg_kwargs)
        return {p.element: p.symbol for p in pmgset.potcar}

    def get_potcar_family(self) -> str:
        """
        Get the POTCAR family used by the input set.

        Retrieves the pseudopotential functional family from the pymatgen
        input set configuration. Converts underscore notation to dot notation
        (e.g., PBE_54 becomes PBE.54) for aiida-vasp compatibility.

        :returns: Name of the POTCAR family
        :rtype: str
        """
        return self._pmg_class.CONFIG['POTCAR_FUNCTIONAL'].replace('_', '.')

    def get_kpoints(self, structure: orm.StructureData) -> Optional[orm.KpointsData]:
        """
        Return a KpointsData object for the given structure.

        Converts the k-point specification from the pymatgen input set to an
        aiida-vasp compatible KpointsData object. Supports Gamma-centered,
        Monkhorst-Pack, and automatic k-point generation modes.

        :param structure: Crystal structure for k-point generation
        :type structure: orm.StructureData

        :returns: K-points data object, or None if no k-points are specified
        :rtype: orm.KpointsData or None
        """
        ps = structure.get_pymatgen()
        pmgset = self._pmg_class(ps, **self._pmg_kwargs)
        if pmgset.kpoints is None:
            return None
        # Currently only supports Gamma and Monkhorst-Pack
        kpoints_data = pmg_kpoints2kpointsdata(pmgset.kpoints, structure)
        return kpoints_data

    def get_kpoints_spacing(self, structure: orm.StructureData) -> Optional[float]:
        """
        Get the k-point spacing used by the input set.

        Extracts the KSPACING parameter from the pymatgen input set and converts
        it to the format expected by aiida-vasp (dividing by 2π).

        :param structure: Crystal structure for the calculation
        :type structure: orm.StructureData

        :returns: K-point spacing value or None if not specified
        :rtype: float or None
        """
        ps = structure.get_pymatgen()
        pmgset = self._pmg_class(ps, **self._pmg_kwargs)
        incar_dict = {key.lower(): value for key, value in pmgset.incar.items()}
        kspacing = incar_dict.pop('kspacing', None)
        if kspacing is not None:
            return kspacing / np.pi / 2
        return None


def pmg_kpoints2kpointsdata(pmg_kpoints, structure: orm.StructureData) -> orm.KpointsData:
    """
    Convert a pymatgen Kpoints object to an AiiDA KpointsData object.

    Handles conversion between different k-point generation modes:
    - Gamma-centered grids
    - Monkhorst-Pack grids (with appropriate shift corrections)
    - Automatic k-point generation based on spacing

    :param pmg_kpoints: Pymatgen Kpoints object to convert
    :type pmg_kpoints: pymatgen.io.vasp.inputs.Kpoints
    :param structure: AiiDA structure data for setting the unit cell
    :type structure: orm.StructureData

    :returns: Converted k-points data object
    :rtype: orm.KpointsData

    :raises ValueError: If the k-point style is not supported
    """

    # Currently only supports Gamma and Monkhorst-Pack
    style = pmg_kpoints.style
    mesh = pmg_kpoints.kpts[0]
    kpoints = orm.KpointsData()

    if style == KpointsSupportedModes.Gamma:
        # aiida-vasp defaults to use Gamma-centering mode when constructing the KPOINTS file
        gamma_shifts = (0, 0, 0)
    elif style == KpointsSupportedModes.Monkhorst:
        # If a MP grid is supplied, we just add a -0.5 shift to the gamma-centering grid to make it
        # equivalent to a MP centred grid
        # See https://www.vasp.at/wiki/index.php/KPOINTS
        gamma_shifts = []
        for i in mesh:
            if i % 2 == 1:
                # Odd division - do nothing
                gamma_shifts.append(0)
            else:
                # Even division - add shifts
                gamma_shifts.append(-0.5)
    elif style == KpointsSupportedModes.Automatic:
        # The automatic mode is used with a length R_k
        # We convert it to a mesh based on the structure's lattice
        # See: https://www.vasp.at/wiki/index.php/KPOINTS#Automatic_k-point_mesh
        kspacing = 2 * np.pi / pmg_kpoints.kpts[0][0]
        kpoints.set_cell_from_structure(structure)
        kpoints.set_kpoints_mesh_from_density(kspacing)
        return kpoints
    else:
        raise ValueError(f'Unsupported kpoint style: {style}')
    # Using explicit meshes
    shifts = pmg_kpoints.kpts_shift
    # Construct AiiDA KpointsData object
    kpoints.set_kpoints_mesh(mesh, offset=[i + j for i, j in zip(shifts, gamma_shifts)])
    kpoints.set_cell_from_structure(structure)
    return kpoints
