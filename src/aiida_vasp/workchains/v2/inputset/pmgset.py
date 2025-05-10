"""
Module for using pymatgen.io.vasp.sets based input sets.
"""

from typing import Union

import numpy as np
from aiida import orm

from .base import InputSet


class PymatgenInputSet(InputSet):
    """
    Input set using pymatgen.io.vasp.sets.
    Provides basic compatibility with pymatgen sets.
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

    def __init__(self, set_name: str, overrides=None, verbose=None, pmg_kwargs=None):
        """
        Instantiate a PymatgenInputSet.

        :param overrides: A dictionary of overrides for the input set.
        :param verbose: If True, print additional information.
        :param pmg_kwargs: A dictionary of keyword arguments to pass to the pymatgen input set.
        """
        assert set_name in self.KNOWN_SETS, f'Unsupported set name: {set_name}'
        super().__init__(set_name, overrides=overrides, verbose=verbose)
        self._pmg_kwargs = pmg_kwargs or {}

    def _load_data(self):
        """Load the data for the input set."""
        try:
            import pymatgen.io.vasp.sets as pmg_sets
        except ImportError:
            raise ImportError('pymatgen has to be installed to use this feature.')

        self._pmg_class = getattr(pmg_sets, self.set_name)

    def get_input_dict(self, structure, raw_python=True) -> Union[dict, orm.Dict]:
        """
        Compute the input parameters for a VASP calculation using pymatgen.io.vasp.sets.
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

    def get_pp_mapping(self, structure):
        """Get the potential mapping used by the input set."""
        ps = structure.get_pymatgen()
        pmgset = self._pmg_class(ps, **self._pmg_kwargs)
        return {p.element: p.symbol for p in pmgset.potcar}

    def get_potcar_family(self):
        """
        Get the POTCAR family used by the input set.
        Note that we replace `_` by `.`, so PBE_54 is PBE.54
        """
        return self._pmg_class.CONFIG['POTCAR_FUNCTIONAL'].replace('_', '.')

    def get_kpoints(self, structure) -> orm.KpointsData:
        """
        Return a KpointsData object for the given structure.
        """
        ps = structure.get_pymatgen()
        pmgset = self._pmg_class(ps, **self._pmg_kwargs)
        if pmgset.kpoints is None:
            return None
        # Currently only supports Gamma and Monkhorst-Pack
        kpoints_data = pmg_kpoints2kpointsdata(pmgset.kpoints)
        kpoints_data.set_cell_from_structure(structure)
        return kpoints_data

    def get_kpoints_spacing(self, structure):
        """Get the spacing of the kpoints used by the input set."""
        ps = structure.get_pymatgen()
        pmgset = self._pmg_class(ps, **self._pmg_kwargs)
        incar_dict = {key.lower(): value for key, value in pmgset.incar.items()}
        kspacing = incar_dict.pop('kspacing', None)
        if kspacing is not None:
            return kspacing / np.pi / 2
        return None


def pmg_kpoints2kpointsdata(pmg_kpoints):
    """Convert a pymatgen Kpoints object to an AiiDA KpointsData object."""
    from pymatgen.io.vasp.inputs import KpointsSupportedModes

    # Currently only supports Gamma and Monkhorst-Pack
    style = pmg_kpoints.style
    mesh = pmg_kpoints.kpts[0]
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
    else:
        raise ValueError(f'Unsupported kpoint style: {style}')
    shifts = pmg_kpoints.kpts_shift
    # Construct AiiDA KpointsData object
    kpoints = orm.KpointsData()
    kpoints.set_kpoints_mesh(mesh, offset=[i + j for i, j in zip(shifts, gamma_shifts)])
    return kpoints
