---
file_format: mystnb
kernelspec:
  display_name: Python 3
  name: python3
---
(silicon_sp_tutorial)=


# Single point and general calculation


:::{note}
This tutorial is can be executed as a jupyter notebook.
Please copy the `mock_registry` folder found to your directory of execution, otherwise the mock
code cannot locate the pre-computed calculations to be used as dummy output.

This notebook can be downloaded as **{nb-download}`silicon_sp.ipynb`** and {download}`silicon_sp.md`
:::

`VaspBuilderUpdater` provides a simplified interface for setting up calculations
using pre-define default inputs for calculation and for workflow execution.

In jupyter notebook, we normally need to load the necessary AiiDA environment:


```
%load_ext aiida
%aiida
```

However, for this tutorial, we will use a temporary AiiDA profile.
All results will be destroyed once the session close.

We also assume that you have already installed the aiida-vasp plugin and the VASP has been installed locally on your computer.

If you do not have VASP available, run the code cell below to activate `mock-vasp` which will
use cached results instead of running actual VASP calculations.

## Setting up the basic environment

Before we start any calculation, we need to do a few basic setups:
1. Setup a `Computer` node, which is this computer (localhost)
2. Tell AiiDA where to find the VASP executable, or using the `mock-vasp` executable.
3. Upload the pseduopotentials (POTCAR) family

The following code creates  `Computer` node:

```{code-cell} python3
from aiida_vasp.utils.temp_profile import *
print(load_temp_profile())


# Uncomment the below line to create a localhost Computer if you have not done so
comp = orm.Computer('localhost', 'localhost', transport_type='core.local', scheduler_type='core.direct')
comp.store()


# Some configuration may be needed for first-time user
comp.set_workdir('/tmp/aiida_run/')
comp.configure()
```

Here we use `mock-vasp` to simulate the VASP executable.

```{code-cell} python3
from pathlib import Path
import os
vasp_path = !which mock-vasp
vasp_code = orm.InstalledCode(comp, vasp_path[0], default_calc_job_plugin='vasp.vasp')
print(vasp_path[0])
vasp_code.label ='mock-vasp'
vasp_code.store()
os.environ['MOCK_VASP_REG_BASE'] = str((Path() / 'mock_registry').absolute())
os.environ['MOCK_VASP_UPLOAD_PREFIX'] = 'singlepoint'
print(os.environ['MOCK_VASP_REG_BASE'])
```
If you have VASP installed, uncomment and run the the code below to create the `InstalledCode` node`:

```{code-cell} python3
#vasp_path = !which vasp_std
#vasp_code = orm.InstalledCode(comp, vasp_path[0], default_calc_job_plugin='vasp.vasp')
#vasp_code.label ='vasp'
#vasp_code.store()
```

:::{hint}
Setting the `MOCK_VASP_VASP_CMD` environment variable will allow the mock-vasp to run real VASP
calculation and add the results to the registry whenever needed.
This environmental variable should be set when generating the test/demo data.
:::

We also need to upload the pseudopotential family, here we use the sample Si POTCAR with a family name of `PBE.EXAMPLE`:

```{code-cell} python3
from aiida_vasp.data.potcar import PotcarData, PotcarFileData
from pathlib import Path

print(PotcarData.upload_potcar_family(str(Path('potcars').absolute()), "PBE.EXAMPLE", "PBE.EXAMPLE"))
```


## Running the calculation

We can now set up the calculation, but first we need to create a `StructureData` node for the Si structure. Here we use the `ase` package to create the structure.

```{code-cell} python3
from ase.build import bulk
si = bulk('Si', 'diamond', 5.4)
si_node = orm.StructureData(ase=si)
```

:::{tip}
It is also possible to create a `StructureData` node from a `pymatgen.core.Structure` object.

```python
from pymatgen.core import Structure
si = Structure.from_spacegroup('Fm-3m', Lattice.cubic(5.4), ['Si'], [[0, 0, 0]]).get_primitive_structure()
si_node = orm.StructureData(pymatgen=si)
```
:::


using the `VaspBuilderUpdater` class:

```{code-cell} python3
from aiida import orm
from aiida_vasp.workchains.v2 import VaspBuilderUpdater

# This instantiate a VaspBuilderUpdater object and apply the preset
# The default name is VaspPreset stored in the code repository.
# You can place your own preset at ~/.aiida-vasp/ and use them for production
# calculations.
upd = VaspBuilderUpdater().apply_preset(si_node, code='mock-vasp@localhost')
# This preset defaults to use the PBE.54 family, we switch to 'PBE.EXAMPLE' family
# just uploaded.
upd.builder.potential_family = 'PBE.EXAMPLE'
```

The code block above create a `VaspBuilderUpdater` object and apply the preset for the Si structure.
We can verify that this configures the `ProcessBuilder` object with the correct inputs:

```{code-cell} python3
upd.builder
```

If this looks all fine, we can run the calculation using the `run_get_node` method:

```{code-cell} python3
results = upd.run_get_node()
results
```

## Accessing the results

The returned `results` contains the outputs as well as the `WorkChainNode` object representing the workflow that was executed.

```{code-cell} python3
workchain_node = results.node
workchain_node
```

Note that our `workchain_node` has a *pk* as well as a *uuid*. Both of them are identifiers for the node. You can load the node using the `load_node` method:

```python
from aiida.orm import load_node
node = load_node(workchain_node.pk)
```

at a later time to access the results.

:::{hint}
Although *pk* and *uuid* are both unique identifiers for this node, they serve different purposes. The *pk* is a unique identifier within the *current aiida database*, while the *uuid* is a truelly unique identifier that would always refer to the same node, even if the data is exported and imported to different databases.
:::

Useful information such as the total energy, forces and stresses are stored in the `misc` output node:

```{code-cell} python3
workchain_node.outputs.misc.get_dict()
```

The advantage of using a database like AiiDA is that we can easily query for the results of our calculations. For example, to get all calculations that used the `PBE.EXAMPLE` family:

```{code-cell} python3
from aiida.orm import QueryBuilder
from aiida.plugins import WorkflowFactory
q = QueryBuilder()
q.append(orm.WorkChainNode, tag='workchain', project=['*'])
q.append(orm.Str, with_outgoing='workchain', edge_filters={'label': 'potential_family'},
         filters={'attributes.value': 'PBE.EXAMPLE'})
q.all()
```

The `QueryBuilder` object allows us to construct complex queries path that combine different nodes and their properties.
The first `append` method says that we need a `WorkChainNode` and we want to return the node itself.
The second `append` method defines that such `WorkChainNode` should have an outgoing `Dict` node with a `label` of `potential_family` and a `value` of `PBE.EXAMPLE`. This filters the results to only include nodes that used the `PBE.EXAMPLE` family.
The final `q.all()` method returns a list of all the combinations of node and edges that match the query.
