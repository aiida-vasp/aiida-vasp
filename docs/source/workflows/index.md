---
file_format: mystnb
kernelspec:
  display_name: Python 3
  name: python3
substitutions:
  VaspWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.vasp.VaspWorkChain>`"
  VaspBandsWorkchain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.bands.VaspBandsWorkchain>`"
  VaspRelaxWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.relax.VaspRelaxWorkChain>`"
  VaspConvergenceWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.converge.VaspConvergenceWorkChain>`"
  calcfunction: "{py:class}`calcfunction <aiida.engine.calcfunction>`"
  workfunction: "{py:class}`calcfunction <aiida.engine.workfunction>`"
---

(using-workflows)=

```{code-cell}
:tags: [hide-cell]

from aiida_vasp.utils.temp_profile import load_temp_profile
load_temp_profile()
```

# Workflows in aiida-vasp

:::{note}
See [this tutorial](#silicon_sp_tutorial) for a quick tour on how to use workflows in aiida-vasp.
:::

There are several workflows bundled with aiida-vasp. They can be referred using the entry point started with `vasp.`

For example, the following code load the standard `VaspWorkChain` in a shell launched by using the command `verdi shell`:

```python
from aiida.plugins import WorkflowFactory  # This can be omitted as it is imported by default with verdi shell
vasp_wc = WorkflowFactory('vasp.v2.vasp')
```

:::{note}
As the continued development of aiida-vasp takes place, the list of workflows may change. The latest workflows stack
are named with the `v2` prefix as they are not compatible with the previous stack.
They will become the default in the next major release of aiida-vasp, after which one can omit the `v2` prefix.
:::

The `VaspWorkChain` is the main workchain that performs a VASP calculation from start to finish.
One can view it as a improved version of of the {{ VaspCalculation }} as it takes care input generation and validation.
It also includes several error handling mechanisms to ensure that the calculation is successful and that the output is valid.
For example, if a geometry optimization run fails to converge due to insufficient wall time requested, the workchain will resubmit a new calculation starting from the last geometry.
The main objective is to ensure the completion of the calculation with the parameters originally specified.

:::{note}
This means that {{ VaspWorkChain }} will not change any parameters that may render the calculated energies incompatible, such as the energy cut off or the k-point grid. However, it may change the electronic solver,
the geometry optimisation algorithm or of the step size.
:::

The {{ VaspWorkChain }} is designed to be general-purpose so it should support any types of VASP calculations.
If it gives *false-positive* assertion of errors, please report them as issues on the [aiida-vasp issue tracker](https://github.com/aiida-vasp/aiida-vasp/issues).
You can also try to turn off the {py:func}`process handler <aiida.engine.process.workchains.restart.valid_handler_overrides>` that raises the error.


## General design principles

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

To print the available settings, one can use:

```{code-cell}
from aiida.plugins import WorkflowFactory
print(WorkflowFactory('vasp.v2.relax').option_class.aiida_description())
```

PS you can also print the input and output ports of the workchain using:

```{code-cell}
from aiida.plugins import WorkflowFactory
!verdi plugin list aiida.workflows vasp.v2.relax
```

By default, every input to the workchain has to be specified in full before submission, this can be quiet tedious for daily calculation.
To simplify the input, we have implemented the `BuilderUpdater` class that can automatically update the builder with default values.

The user may write default values and store them in an YAML file to ensure consistent settings are used across multiple projects.


## Convergence workchain

The {{ VaspConvergenceWorkChain }} is a simple workflow that runs a series of VASP calculations with different parameters and checks if the results converge.
The convergence of cut off energy and kpoints are currently implemented.

As metioned above, the inputs to the {{ VaspWorkChain }} should be placed into the `vasp` namespace.
The convergence settings are specified using the `convergence_settings` input which is a `Dict` containing the following keys:

```
print(WorkflowFactory('vasp.v2.relax').option_class.aiida_description())
```

## Band structure workflow

The {{ VaspBandsWorkchain }} is a workflow for calculating the band structure of a material using VASP.
A band structure typically involves computing the ground state electron density then using this fixed density to
solve for the eigenvalues of the Kohn-Sham equation at specific k-points in the Brillouin zone.

Typically, a path along which the eigenvalues are computed is generated based on the point group symmetry of the
input structure.
There are approaches to generate this path automatically，here we default to using `seekpath`, but it can be
switched to using the paths generated by `sumo`.

Another complication is that the path generated is for a specific primitive-cell configuration (as there are infinite ways of choosing the primitive cell).
Hence, a common mistake is to blindly using the path of the input cell, which may not be the standardized primitive cell.
Here, the workchain handles this internally, and the generated standardized primitive cell is returned by the workchain as one of the outputs.

In addition, an exposed `relax` namespace for running {{ VaspRelaxWorkChain }} exists and the workchain will perform
geometry optimization before the band structure calculation if it is specified.

The parameters for the scf (for generating the charge density) the actual band structure structure calculation should be specified under the exposed {{ VaspWorkChain }} namespace called  `scf` and `bands`.
An additional `dos` namespace is also exposed for calculating the density of states and can be specified if desired.


:::{note}
The `scf` namespace should always be specified, while specifying `bands` namespace is only needed if the
input nodes should be different from that in the `scf` namespace. The same rule applies to the `dos` namespace.
:::

Similar to the {{ VaspRelaxWorkChain }} the behavor of the {{ VaspBandsWorkchain }} can be controlled using the `band_settings` input:

```{code-cell}
from aiida.plugins import WorkflowFactory
print(WorkflowFactory('vasp.v2.bands').option_class.aiida_description())
```


The {{ VaspHybridBandsWorkChain }} is an variant of the {{ VaspBandsWorkchain }} for running band structure calculation with hybrid functional.
In this case, the potential is not completely determined from the electron density, hence one cannot use the standard
approach that first compute the ground state electron density and then use it to solve the Kohn-Sham equation.
Instead, the Kohn-Sham equation has to be solved self-consistently, and the k-points along the path are inserted
as *zero-weighted k-points*.

The {{ VaspHybridBandsWorkchain }} is designed for this purpose.
In addition, the large compute cost of hybrid functional means it may be advantageous to split the full k-point path into smaller sub-paths,
and run multiple self-consistent calculations in parallel instead of doing a single large calculation,
given the constraints of the available computing resources.
The number of kpoints included in each sub-path can be specified using the `kpoints_per_subpath` input.

:::{hint}
Set `kpoints_per_subpath` to a very large number  to run a single self-consistent calculation with all k-points.
:::


[vasp]: https://www.vasp.at
