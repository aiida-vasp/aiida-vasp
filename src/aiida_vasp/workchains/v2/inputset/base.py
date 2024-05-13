"""
Module for preparing standardised input for calculations
"""

from copy import deepcopy
from math import pi
from pathlib import Path

import yaml
from aiida.orm import Dict, KpointsData
from ase import Atoms

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


def get_library_path():
    """Get the path where the YAML files are stored within this package"""
    return Path(__file__).parent


class InputSet:
    """
    Base class representing an inputs set.

    Not useful on its own, should be subclass for convenient definition of inputs
    for high-throughput calculations.
    """

    # path from which the set yaml files are read
    _load_paths = (get_library_path(), Path('~/.inputsets').expanduser())

    def __init__(self, set_name, structure, overrides=None, verbose=True):
        """
        Initialise an InputSet

        Args:
          set_name: Name of the set to be loaded
          structure: Structure used for calculation, can be StructureData or Atoms
          overrides: A dictionary of overriding inputs.

        """
        self.structure = structure
        self.set_name = set_name

        if overrides is None:
            overrides = {}
        self.overrides = overrides

        self._presets = None
        self.verbose = verbose
        self._load_data()

    def get_input_dict(self, raw_python=True):
        """
        Get a input dictionary for VASP
        """

        out_dict = deepcopy(self._presets['global'])

        # Set-per atom properties
        natoms = self.natoms
        for key, value in self._presets.get('per_atom', {}).items():
            out_dict[key] = value * natoms

        self.apply_overrides(out_dict)

        if raw_python:
            return out_dict
        return Dict(dict=out_dict)

    def _load_data(self):
        """Load stored data"""
        set_path = None
        for parent in self._load_paths:
            set_path = parent / (self.set_name + '.yaml')
            if set_path.is_file():
                break
        if set_path is None:
            raise RuntimeError(f'Cannot find input set definition for {self.set_name}')

        if self.verbose:
            print(f'Using input set file at: {set_path}')

        with open(set_path, encoding='utf-8') as fhd:
            self._presets = yaml.load(fhd, Loader=yaml.FullLoader)

    @property
    def natoms(self):
        if isinstance(self.structure, Atoms):
            return len(self.structure)
        return len(self.structure.sites)

    @property
    def elements(self):
        if isinstance(self.structure, Atoms):
            return self.structure.get_chemical_symbols()
        return [kind.name for kind in self.structure.kinds]

    def apply_overrides(self, out_dict):
        """Apply overrides stored in self.overrides to the dictionary passed"""
        for name, value in self.overrides.items():
            # Keys ends with '_mapping' are treated differently here
            # Those valuse should have been applied already implemented in the `get_input_dict` method.
            if '_mapping' in name or '_list' in name:
                continue
            # Delete the key
            if value is None:
                out_dict.pop(name, None)
            else:
                out_dict[name] = value

    def get_kpoints(self, density):
        """
        Return a kpoints object for a given density

        Args:
          density: kpoint density in 2pi Angstrom^-1 (CASTEP convention)

        Returns:
          An KpointsData object with the desired density
        """

        kpoints = KpointsData()
        kpoints.set_cell(self.structure.cell)
        kpoints.set_kpoints_mesh_from_density(density * 2 * pi)
        return kpoints
