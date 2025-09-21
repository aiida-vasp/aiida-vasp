---
file_format: mystnb
kernelspec:
  display_name: Python 3
  name: python3
myst:
  substitutions:
    VaspWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.vasp.VaspWorkChain>`"
    calcfunction: "{py:class}`calcfunction <aiida.engine.calcfunction>`"
    workfunction: "{py:class}`calcfunction <aiida.engine.workfunction>`"
    VaspCalculation: "{py:class}`VaspCalculation <aiida_vasp.calcs.vasp.VaspCalculation>`"
    VaspBandsWorkchain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.bands.VaspBandsWorkchain>`"
    VaspRelaxWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.relax.VaspRelaxWorkChain>`"
    VaspConvergenceWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.converge.VaspConvergenceWorkChain>`"
    hybrid_bands: "{py:class}`VaspHybridBandsWorkChain <aiida_vasp.workchains.v2.bands.VaspHybridBandsWorkchain>`"
    VaspInputGenerator: "{py:class}`VaspInputGenerator<aiida_vasp.protocols.BaseInputGenerator>`"

---
(workflows)=
# Workflows

:::{note}
This notebook can be downloaded as **{nb-download}`workflows.ipynb`** and {download}`workflows.md`
:::

The [Workchain] class is the central piece that enables workflows to be run with aiida-vasp.
By composing one or several [Workchain] classes, one can make a workflow.
As single [WorkChain] class may launch one or several calculations, or it may launch children [WorkChain]s to achieve the designed functionality.

For  any short-running python code, the workchain can run them directly as {{ calcfunction }} or {{ workfunction }} directly, and the provenance will be recorded accordingly.

It is important to note that however, long-running computational tasks *should not* be run directly in the code as it
will delay or block the operation of the [daemon].

We would like to encourage users to build workchains and/or compose existing ones into more advanced workflows that we can all share and benefit from.
You may want to visit [this page](https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/write_workflows.html) to learn more about WorkChains and how to build them.

One should note that the advantage of using a provenance-preserving engine like AiiDA is that you do
not have to define a workflow in order to have the calculations steps recorded and stored.
It is perfectly fine to conduct exploration studies using the basic workchains and use {py:func}`calcfunction <aiida.engine.processes.functions.calcfunction>`  to link the outputs/inputs together for provenance.



## Design principles

The rest of the bundled workchain are designed to run `VaspWorkChain` as the basic unit of work.
This means that they expect error-correction functionalities to be embedded in the `VaspWorkChain` so they
do not need to explicitly handle errors.

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
To simplify the input, we have implemented the {{ VaspInputGenerator }} class that can automatically update the builder with default values.
See [this page](#workflow_inputs) for more information.

The user may write default values and store them in an YAML file to ensure consistent settings are used across multiple projects.


PS you can also print the input and output ports of the workchain using:

```{code-cell}
from aiida.plugins import WorkflowFactory
!verdi plugin list aiida.workflows vasp.v2.relax
```

(bundled_workflows)=
## Workflows included in aiida-vasp


There are several workflows bundled with aiida-vasp. They can be referred using the entry point started with `vasp.`

For example, the following code load the standard `VaspWorkChain` in a shell launched by using the command `verdi shell`:

```python
from aiida.plugins import WorkflowFactory  # This can be omitted as it is imported by default with verdi shell
vasp_wc = WorkflowFactory('vasp.vasp')
```

:::{hint}
You may see something like `vasp.v2.vasp` as entry point in the document:
- The first `vasp` means the entrypoint is from `aiida-vasp` plugin
- The second part `v2` is a version tag, refers to the `v2` version
- The last part `vasp` refers to the `vasp` workchain/calculation included in the plugin.

The latest version of the workchain is selected if the `v2` is omitted. We use this syntax to allow
some backward compatibility during the development.
:::

The {{ VaspWorkChain }} is the main workchain that performs a VASP calculation from start to finish.
One can view it as a improved version of of the {{ VaspCalculation }} as it takes care input generation and validation.
It also includes several error handling mechanisms to ensure that the calculation is successful and that the output is valid.
For example, if a geometry optimization run fails to converge due to insufficient wall time requested, the workchain will resubmit a new calculation starting from the last geometry.
The main objective is to ensure the completion of the calculation with the parameters originally specified.

{{ VaspWorkChain }} will not change any parameters that may render the calculated energies incompatible, such as the energy cut off or the k-point grid. However, it may change the electronic solver,
the geometry optimisation algorithm or of the step size.

The {{ VaspWorkChain }} is designed to be general-purpose so it should support any types of VASP calculations.
If it gives *false-positive* assertion of errors, please report them as issues on the [aiida-vasp issue tracker](https://github.com/aiida-vasp/aiida-vasp/issues).
You can also try to turn off the {py:func}`process handler <aiida.engine.process.workchains.restart.valid_handler_overrides>` that raises the error.



This section we give some brief introduction to the bundled workflows in AiiDA-VASP.

### Convergence workchain

The {{ VaspConvergenceWorkChain }} is a simple workflow that runs a series of VASP calculations with different parameters and checks if the results converge.
The convergence of cut off energy and kpoints are currently implemented.

As metioned above, the inputs to the {{ VaspWorkChain }} should be placed into the `vasp` namespace.
The convergence settings are specified using the `convergence_settings` input which is a `Dict` containing the following keys:

```python
print(WorkflowFactory('vasp.v2.relax').option_class.aiida_description())
```

### Relaxation workchain

The {{ VaspRelaxWorkChain }} is a simple workflow that runs a VASP relaxation calculation.
It will run VASP geometry optimizations until the specified convergence criteria are met.

This may involve one or more actual VASP calculations. This is because:

- A single VASP calculation may not fully relax the structure, especially when the maximum number of ionic steps is set to a relatively small value.
- For variable cell geometry optimization, multiple VASP calculations are required as each restart resets the basis set, otherwise the effective cut off energy can change.
- A final singlepoint calculation may be needed to ensure that the energy is consistent with the cut off, if the lattice has been changed.

The inputs to the {{ VaspRelaxWorkChain }} should be placed into the `vasp` namespace.
The convergence settings are specified using the `relax_settings` input which is a `Dict` containing the following keys:

```{code-cell}
from aiida.plugins import WorkflowFactory
print(WorkflowFactory('vasp.v2.relax').option_class.aiida_description())
```

Note the keys such as `algo`, `steps`, `force_cutoff` are translated into INCAR tags (`IBRION`, `NSW`, `EDIFFG`, etc.), so one should not explicitly set these tags in the `parameters` input.

:::{hint}
This means one can quickly reuse the `parameters` from a single point calculation for a relaxation and *vice versa*.
:::

See [this tutorial](#silicon_relax) for an example of how to run the {{ VaspRelaxWorkChain }}.

### Band structure workflow

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


The {{ hybrid_bands }} is an variant of the {{ VaspBandsWorkchain }} for running band structure calculation with hybrid functional.
In this case, the potential is not completely determined from the electron density, hence one cannot use the standard
approach that first compute the ground state electron density and then use it to solve the Kohn-Sham equation.
Instead, the Kohn-Sham equation has to be solved self-consistently, and the k-points along the path are inserted
as *zero-weighted k-points*.

The {{ hybrid_bands }} is designed for this purpose.
In addition, the large compute cost of hybrid functional means it may be advantageous to split the full k-point path into smaller sub-paths,
and run multiple self-consistent calculations in parallel instead of doing a single large calculation,
given the constraints of the available computing resources.
The number of kpoints included in each sub-path can be specified using the `kpoints_per_subpath` input.

:::{hint}
Set `kpoints_per_subpath` to a very large number  to run a single self-consistent calculation with all k-points.
:::

See [this tutorial](#band_dos) for an example of how to run the {{ VaspBandsWorkchain }}.

[vasp]: https://www.vasp.at
[workchain]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/workflows.html#work-chains
[daemon]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/topics/daemon.html
