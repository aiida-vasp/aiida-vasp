---
file_format: mystnb
kernelspec:
  display_name: Python 3
  name: python3
execution:
  timeout: 300
---
(silicon_relax)=

# Geometry optimisation

:::{note}
This notebook can be downloaded as **{nb-download}`silicon_relax.ipynb`** and {download}`silicon_relax.md`
:::

In this example we will perform a geometry optimisation of silicon using VASP.
We will use the VASP code and the AiiDA plugin for VASP.

## Setting up  the environment

The code block below configures the environment for this example.
Please see the [single point example](#silicon_sp_tutorial) for more details.


```{code-cell} python3
from aiida_vasp.utils.temp_profile import *
print(load_temp_profile())


# Uncomment the below line to create a localhost Computer if you have not done so
comp = orm.Computer('localhost', 'localhost', transport_type='core.local', scheduler_type='core.direct')
comp.store()


# Some configuration may be needed for first-time user
from pathlib import Path
import os
comp.set_workdir('/tmp/aiida_run/')
comp.configure()
vasp_path = !which mock-vasp
vasp_code = orm.InstalledCode(comp, vasp_path[0], default_calc_job_plugin='vasp.vasp')
vasp_code.label ='mock-vasp'
vasp_code.store()
os.environ['MOCK_VASP_REG_BASE'] = str((Path() / 'mock_registry').absolute())
os.environ['MOCK_VASP_UPLOAD_PREFIX'] = 'relax'

# Upload the POTCAR files
from aiida_vasp.data.potcar import PotcarData, PotcarFileData
from pathlib import Path

print(PotcarData.upload_potcar_family(str(Path('potcars').absolute()), "PBE.EXAMPLE", "PBE.EXAMPLE"))

# Setting up the silicon structure
from ase.build import bulk
si = bulk('Si', 'diamond', 5.4)
si_node = orm.StructureData(ase=si)
```

## Running the relaxation

Similar to the single point calculation tutorial, we will use a `BuilderUpdater` to setup the inputs
for the `VaspRelaxWorkChain`.

```{code-cell}
from aiida import orm
from aiida_vasp.workchains.v2 import VaspRelaxUpdater

upd = VaspRelaxUpdater().apply_preset(si_node, code='mock-vasp@localhost')
upd.builder.vasp.potential_family = 'PBE.EXAMPLE'
upd.builder
```

Here we can see that the input nodes are different from that of the single point calculation (workflow).
Nodes such as the `settings` and `parameters` go into a `vasp` input names rather than at the root level.
This is because the `VaspRelaxWorkChain` exposes the inputs of the `VaspWorkChain` in the `vasp` namespace.

Also note that the `structure` input is still at the `root` level.
This is because the atomic structure is an essential input for the relaxation (and another calculations).

Note that there is a also a `relax_settings` input, which contains control parameters for the relaxation.
You can see the available options by inspecting the field of the builder:

```{code-cell}
upd.get_input_help('relax_settings')
```

or simply hit `<Shift>+<Tab>` after `upd.builder.relax_settings` in the Notebook code cell.

We can now run the relaxation using the same `run_get_node` method as in the single point example.

```{code-cell}
:tags: [remove-stderr]
results = upd.run_get_node()
```

The `VaspRelaxWorkChain` also has a `misc` output node, which is in fact the output of the last `VaspWorkChain` it launched.

```{code-cell}
results.node.outputs.misc.get_dict()
```

The final relaxed structure is stored in the `relax.structure` output node.

```{code-cell}
relaxed_node = results.node.outputs.relax.structure
print(f"Volume before relaxations: {si_node.get_cell_volume():3f} A^3")
print(f"Volume after relaxations: {relaxed_node.get_cell_volume():3f} A^3")
```

## Why use `VaspRelaxWorkChain`?

Surely VASP has the option to run relaxation itself (`ISIF=3`, `IBRION=2`, `NSW=60`), so why do we need a `VaspRelaxWorkChain` in the first place?
The reason is tha `VaspRelaxWorkChain` run more strict check and verifies if the structure is really fully relaxed.
Also, since VASP run variable cell relaxation in the *fixed-basis* mode, the effective plane wave cut off energies changes with the cell volume inside a single calculation.
This means multiple calculations are needed to get fully (and self-consistently) relaxed structures.
In addition, the `VaspRelaxWorkChain` will run a single point calculation in the end to get the correct energy that is consistent with the plane wave cut-off energy (`ENCUT`) originally set.

The `called` and `called_descendants` property can be used to inspect the sub-processes launched by the `VaspRelaxWorkChain`

```{code-cell}
print("Directly called workchains:", results.node.called)
print("All called processes:", results.node.called_descendants)
```
