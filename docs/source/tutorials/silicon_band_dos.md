---
file_format: mystnb
kernelspec:
  display_name: Python 3
  name: python3
execution:
  timeout: 300
---
(band_dos)=

# Band structure and density of states

:::{note}
This notebook can be downloaded as **{nb-download}`silicon_band_dos.ipynb`** and {download}`silicon_band_dos.md`
:::

In this example we will perform band structure and DOS calculation for silicon using VASP. We will use the `VaspBandWorkChain` from the `aiida-vasp` plugin.

It is recommended to go through the [single point calculation tutorial](#silicon_sp) first before proceeding with this example.

## Setting up the environment


The code block below configures the environment for this example.
Please see the [single point example](#silicon_sp_tutorial) for more details.


```{code-cell} python3
from aiida_vasp.utils.temp_profile import *
from subprocess import check_output
print(load_temp_profile())


# Uncomment the below line to create a localhost Computer if you have not done so
comp = orm.Computer('localhost', 'localhost', transport_type='core.local', scheduler_type='core.direct')
comp.store()


# Some configuration may be needed for first-time user
import os
from pathlib import Path
comp.set_workdir('/tmp/aiida_run/')
comp.configure()
vasp_path = check_output(['which', 'mock-vasp'], universal_newlines=True).strip()

vasp_code = orm.InstalledCode(comp, vasp_path, default_calc_job_plugin='vasp.vasp')
vasp_code.label ='mock-vasp'
vasp_code.store()
os.environ['MOCK_VASP_REG_BASE'] = str((Path() / 'mock_registry').absolute())
os.environ['MOCK_VASP_UPLOAD_PREFIX'] = 'band'

# Upload the POTCAR files
from aiida_vasp.data.potcar import PotcarData, PotcarFileData
from pathlib import Path

print(PotcarData.upload_potcar_family(str(Path('potcars').absolute()), "PBE.EXAMPLE", "PBE.EXAMPLE"))

# Setting up the silicon structure
from ase.build import bulk
si = bulk('Si', 'diamond', 5.4)
si_node = orm.StructureData(ase=si)
```

## Setting up the band structure and DOS calculation

Similar to the single point calculation tutorial, we will use a `BuilderUpdater` to setup the inputs
for the `VaspBandsWorkChain`:

```{code-cell}
from aiida_vasp.workchains.v2.bands import BandOptions
from aiida_vasp.workchains import VaspBandUpdater

upd = VaspBandUpdater().apply_preset(si_node, code='mock-vasp@localhost')
upd.builder.scf.potential_family = 'PBE.EXAMPLE'
```

The workchain can be modified with several options. These options are stored in the the `band_settings` input node which of the type `orm.Dict`.
The available options can be printed using the `aiida_description()` method.

```{code-cell}
opt = BandOptions()
print(opt.aiida_description())
```

:::{hint}
This can also be done for other Option classes such as `relax_settings` for the `VaspRelaxWorkChain` which uses the `RelaxOptions` class.
:::

## Run and inspect the results

We can now run the workchain and get the returned `WorkChainNode` object.

```{code-cell}
:tags: [remove-stderr]
band_out = upd.run_get_node().node
```

The computed band structure is stored as a `BandsData` object in the `band_structure` output port. We can plot the band structure using the `show_mpl` method provided.

```{code-cell}
band_out.outputs.band_structure.show_mpl()
```
