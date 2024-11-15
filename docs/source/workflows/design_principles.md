---
file_format: mystnb
kernelspec:
  display_name: Python 3
  name: python3
myst:
  substitutions:
    VaspWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.vasp.VaspWorkChain>`"
    VaspBandsWorkchain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.bands.VaspBandsWorkchain>`"
    VaspRelaxWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.relax.VaspRelaxWorkChain>`"
    VaspConvergenceWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.converge.VaspConvergenceWorkChain>`"
    calcfunction: "{py:class}`calcfunction <aiida.engine.calcfunction>`"
    workfunction: "{py:class}`calcfunction <aiida.engine.workfunction>`"
---

```{code-cell}
:tags: [hide-cell]

from aiida_vasp.utils.temp_profile import load_temp_profile
load_temp_profile()
```


# Design principles

The rest of the bundled workchain are designed to run `VaspWorkChain` as the basic unit of work.
This means that they expect error-correction functionalities to be embedded in the `VaspWorkChain` so they
doe not need to explicitly handle errors.

We use the `expose_input` and `expose_outputs` methods of the `WorkChain` class to expose the inputs and outputs of the `VaspWorkChain`.

For example, the inputs to the relax workchain looks like this:

```
VaspRelaxWorkChain
|
|- structure (StructureData of the input structure)
|- vasp (exposed VaspWorkChain inputs)
|- static_calc_settings (settings to override for the static calculation)
|- static_calc_options (options to override for the static calculation)
|- static_calc_parameters (parameters to override for the static calculation)
|- relax_settings (settings controlling the relaxation)
|- verbose
```

Where the inputs specific to the {{ VaspWorkChain }} to be launched as nested inside the `vasp` namespace.
For example, to set the parameters one can use do the following:

```python
from aiida.plugins import WorkflowFactory
builder = WorkflowFactory('vasp.v2.relax').get_builder()
builder.vasp.parameters = Dict(dict={'incar': {'encut': 500, 'isif': 2, 'nsw': 5, 'potim': 0.01}})
```

while when using {{ VaspWorkChain }} directly, one can use:

```python
from aiida.plugins import WorkflowFactory
builder = WorkflowFactory('vasp.v2.vasp').get_builder()
builder.parameters = {'incar': {'encut': 500, 'isif': 2, 'nsw': 5, 'potim': 0.01}}  # This gets converted to a Dict automatically
```

The other options at the top level are specific to the workchain and are used to control its behavior.

The `relax_settings` input is a `Dict` that contains the settings for the relaxation.
These settings are validated at the submission time using the `pydantic` library.

To see the available settings, one can use:

```{code-cell}
from aiida.plugins import WorkflowFactory
opt = WorkflowFactory('vasp.v2.relax').option_class
# opt.<tab> to see all available options
print(opt.aiida_description())
```

By default, every input to the workchain has to be specified in full before submission, this can be quiet tedious for daily calculation.
To simplify the input, we have implemented the [`BuilderUpdater`] class that can automatically update the builder with default values.
See [this page](#workflow_inputs) for more information.

The user may write default values and store them in an YAML file to ensure consistent settings are used across multiple projects.


PS you can also print the input and output ports of the workchain using:

```{code-cell}
from aiida.plugins import WorkflowFactory
!verdi plugin list aiida.workflows vasp.v2.relax
```
