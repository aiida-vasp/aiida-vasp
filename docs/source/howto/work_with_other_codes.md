---
file_format: mystnb
kernelspec:
  name: python3
  display_name: 'Python 3'
myst:
  substitutions:
    get_sumo_bands_plotter: "{py:func}`get_sumo_bands_plotter <aiida_vasp.utils.sumo.get_sumo_bands_plotter>`"
    get_sumo_dos_plotter: "{py:func}`get_sumo_dos_plotter<aiida_vasp.utils.sumo.get_sumo_dos_plotter>`"
---

(other_codes)=
# How to work with other codes

:::{note}
This notebook can be downloaded as **{nb-download}`work_with_other_codes.ipynb`** and {download}`work_with_other_codes.md`
:::

In this guide, we will learn how to work with other code using aiida-vasp.

## How to work with Atomic Simulation Environment (ASE)

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


## How to work with Pymatgen


The `pymatgen.core.Structure` can be converted to `aiida.orm.StructureData` and vice versa.

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

If you have a pymatgen installation with VASP POTCARs configured (following this [guide](https://pymatgen.org/installation.html#potcar-setup)), you can use the `aiida-vasp potcar upload-from-pymatgen` command to upload them to the AiiDA database.
As in the normal upload, the family name must be specified.
Pymatgen distinguishes different POTCAR sets as different *functionals*, so the functional must also be specified.
The `PBE.54` family mentioned in the documentation refers to the `potpaw.54` POTCAR set, which is the `PBE_54` functional as in pymatgen.

:::{note}
Pymatgen defaults to the `PBE` POTCAR set (*functional*) which is quite OLD had been superseded by multiple updated sets .
Certain POTCARs in this set can be problematic (such as the `W_pv`, which is removed in `PBE_54`).
One should avoid using this set unless direct comparison of raw energies with the Materials Project is required.
:::



(pymatgen-vasp-io)=
### Coming from pymatgen based workflows

If you are coming from using pymatgen for setting up VASP input files,
the {py:class}`VaspInputGenerator<aiida_vasp.protocols.BaseInputGenerator>` interface would feel
very familiar,
which uses the same approach of using a pre-defined set of input parameters.

Consider the following code using pymatgen to set up a VASP calculation:

```{code-cell}
:tags: [remove-stderr]
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.util.testing import PymatgenTest

incar_dict = { 'EDIFFG': -1e-2, 'IVDW': 11, 'ISYM':2,'NSW':1500, 'ENCUT':520}
# Load structure from some file
#structure = Structure.from_file("Al_empty.cif")
# We use a built-in structure for now
structure = PymatgenTest.get_structure("CsCl")
inputset = MPRelaxSet(structure = structure,user_incar_settings=incar_dict,
                       user_kpoints_settings={'length':25})
inputset.write_input(output_dir='./DFT_calc',include_cif=True)
inputset
```

Which loads the `Al_empty.cif` file, sets up a `MPRelaxSet` with some user defined settings, and writes the input files to the `DFT_calc` directory.

To achieve a similar (but not equivalent) effect with aiida-vasp:

```{code-cell}
:tags: [remove-cell]
# Configure the temp profile environment
from aiida_vasp.utils.temp_profile import load_temp_profile_with_mock
from aiida_vasp.data.potcar import PotcarData
from pathlib import Path
load_temp_profile_with_mock()
PotcarData.upload_potcar_family(str(Path('../tutorials/potcars').absolute()), "PBE.EXAMPLE", "PBE.EXAMPLE")
```

```{code-cell}
:tags: [remove-stderr]
from aiida import orm
from pymatgen.util.testing import PymatgenTest
from pymatgen.core import Structure
from aiida_vasp.protocols.generator import VaspInputGenerator

#structure = Structure.from_file("Al_empty.cif")
structure = PymatgenTest.get_structure("CsCl")

upd = VaspInputGenerator(protocol="MPRelaxSet")
upd.get_builder(structure=orm.StructureData(pymatgen=structure), code='mock-vasp@localhost', overrides={"potential_family": "PBE.EXAMPLE"})
upd.set_resources(num_machines=1, tot_num_mpiprocs=16)
upd.set_options(max_wallclock_seconds=3600)
upd.builder
```

There are a few differences to note:

1. The `VaspInputGenerator` class is used to set up the input parameters with presets, here we used the `MPRelaxSet` as the set name. The code will automatically load use the `pymatgen.io.vasp.sets.MPRelaxSet` class to setup the input parameters subsequently.
2. In addition to the calculation input, one needs to define resources requested from the computing cluster's scheduler. This is because the `submit` method submits all calculation data to the daemon which takes care the rests, rather than having the user manually transfer the data to the remote machine, submit the job, and then retrieve the results. In fact, what gets submitted is a *workflow* which may apply automatic restarts and error corrections if needed.
3. Care should be taken to valid the **actual** calculation parameters as `MPRelaxSet` returns some parameters that are controllbed by higher-level workchain in the framework of `aiida-vasp`, such as `ibrion`, `nsw` and `isif`. These parameters may need to be removed (set to `None`) via overrides.

The `VaspInputGenerator` also takes an argument of the **preset** name which gives a higher level of control over how the calculation
should be configured. The **preset** includes which [protocol](../concepts/protocols)  should be used, what overrides should be applied as well as how they should be adapted for different types of workflow as well as for different Code/Computers.
For example, different `NCORE` may be applied when running VASP on different machines.
Users may want to define their own **preset** rather than creating/modifying the input sets directly.

Since AiiDA can store the structure files as `StructureData` nodes, it is possible to first read files and then use a single structure as inputs to multiple subsequent calculations, rather than creating new but identical ` StructureData` nodes each time that a calculation is submitted.

:::{note}
We recommend using `VaspRelaxGenerator` rather than `VaspInputGenerator` for running geometry optimisations, as the former offers more checks and control for process.
:::

## How to work with `sumo`

[Sumo](https://smtg-bham.github.io/sumo/) is a code for plotting publication quality electronic band structures and density-of-states.
It can be used with calculations done by aiida-vasp. There are two ways to use sumo with this package:

1. We can export the calculation with `aiida-vasp tools export <node> <folder>` and then use sumo's command line interface to plot.
   This approach works best for DOS plots and for band structure calculations the exported KPOINTS files currently does not have band labels.

2. Use the `vasp.utils.sumo` module to plot the band structure from a band structure workflow with the get_sumo_dos_plotter and get_sumo_bands_plotter functions:
   ```python
   from aiida_vasp.utils.sumo import get_sumo_bands_plotter

   band_node = load_node(<uuid_pk_of_a_bandstructure_calculation>)
   band_data = band_node.outputs.band_structure
   structure = band_node.outputs.primitive_structure
   plotter = get_sumo_dos_plotter(band_data)
   plotter.get_plot()
   ```
   This approach is more convenient when working in a Jupyter Lab environment.

There are other functions such as {py:func}`get_pmg_bandstructure <aiida_vasp.utils.sumo:get_pmg_bandstructuree`, {py:func}`get_sumo_dos_plotter <aiida_vasp.utils.sumo:get_sumo_dos_plotter>` and {py:func}`bandstats <aiida_vasp.utils.sumo:bandstats>` which may be useful when analyzing electronic structures.
