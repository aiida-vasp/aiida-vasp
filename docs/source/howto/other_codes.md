(other_codes)=
# Working with other codes

This section will show how to use aiida-vasp to work with other codes.

## Atomic Simulation Environment (ASE)

ASE can be used to setup input structure or analyse the outputs. The AiiDA core package
already contains routines for converting `ase.Atoms` between `aiida.orm.StructureData`

```python
from aiida import orm
structure = orm.StructureData(ase=ase_atoms)
atoms = structure.get_ase()
```

In addition, ASE can be very useful for visualizing the structures in a notebook environment.

```python
from ase.visualize import view
view(structure_node.get_ase())
```

:::{hint}
You may want to construct a `view` function in your ipython start file (`~/.ipython/profile_default/startup/ase_view.py`) to make it easier to visualize structures with automatic conversion.

```python
try:
    import pymatgen
    import ase
except ImportError:
    pass
else:
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.core import Structure
    from ase.visualize import view as aview

    def view(atoms, *args, **kwargs):
        """
        Allow viewing pymatgen.core.Structure using ase
        """
        if isinstance(atoms, (list, tuple)):
            if isinstance(atoms[0], Structure):
                atoms = [AseAtomsAdaptor.get_atoms(s) for s in atoms]
        elif isinstance(atoms, Structure):
            atoms = AseAtomsAdaptor.get_atoms(atoms)
        elif hasattr(atoms, 'get_ase'):   # check if a aiida.orm.StructureData
            atoms = atoms.get_ase()
        return aview(atoms, *args, **kwargs)
```

:::


## Pymatgen

### Structure conversion

Similarly, pymatgen can be used to convert `pymatgen.core.Structure` to `aiida.orm.StructureData` and vice versa.

```python
from aiida import orm
structure_node = orm.StructureData(pymatgen=pymatgen_structure)
pymatgen_structure = structure_node.get_pymatgen()
```

Then pymatgen's analysis and visualization tools can be used as usual.

### VASP IO with pymatgen

Pymatgen also has its own classes for working with VASP calculations. These objects are not directly supported as inputs to VASP calculation in AiiDA-VASP.
However, it is possible to load these objects from a finished `VaspCalculation` or `VaspWorkChain` and use them in subsequent analysis.

```
from aiida_vasp.utils.pmg import PymatgenAdapter

vasp_calc = load_node('<uuid>')
adapt = PymatgenAdapter(vasp_calc)
vasprun = adapt.vasprun  # Retrieve the pymatgen Vasprun object
```

This is possible by AiiDA-VASP preserves the original calculation output files in the storage.
Behind the scene, the calculation folder is reconstructed inside a temporary directory.

Since exporting the raw files can be slow, cache has been implemented so it is possible to get objects without re-exporting the files every time a property is accessed.
The caches stores the output of the `as_dict` of the corresponding python objects as the `extras` of the calculation node.
Some object, cannot be reconstructed due to the limitations in pymatgen, but they can still be accessed as dictionaries with the property name suffixed with `_dict`.

```
from aiida_vasp.utils.pmg import PymatgenAdapter

vasp_calc = load_node('<uuid>')
# Using with block triggers the cache to be flushed into the storage
with PymatgenAdapter(vasp_calc) as adapt:
    vasprun = adapt.vasprun  # Retrieve the pymatgen Vasprun object

vdict = PymatgenAdapter(vasp_calc).vasprun_dict  # Access the vasprun as a dictionary - this will not export the files again

vasprun PymatgenAdapter(vasp_calc).vasprun  # This WILL re-export the files to the disk and parse using pymatgen again
```

(potcar-from-pymatgen)=
### Uploading pseduopotentials from a pymatgen installation

If you have a pymatgen installation with VASP POTCARs configured, you can use the `verdi data vasp.potcar upload-from-pymatgen` command to upload them to the AiiDA database.
As in the normal upload, the family name must be specified.
Pymatgen distinguishes different POTCAR sets as different *functionals*, so the functional must also be specified.
The `PBE.54` family mentioned in the documentation refers to the `potpaw.54` POTCAR set, which is the `PBE_54` functional as in pymatgen.

:::{note}
Pymatgen defaults to the `PBE` POTCAR set (*functional*) which is quite OLD had been superseded by multiple updated sets .
Certain POTCARs in this set can be problematic (such as the `W_pv`, which is removed in `PBE_54`).
One should avoid using this set unless direct comparison of raw energies with the Materials Project is required.
:::
