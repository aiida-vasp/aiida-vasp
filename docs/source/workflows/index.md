---
myst:
  substitutions:
    VaspWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.vasp.VaspWorkChain>`"
    VaspCalculation: "{py:class}`VaspCalculation <aiida_vasp.calcs.vasp.VaspCalculation>`"
    VaspBandsWorkchain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.bands.VaspBandsWorkchain>`"
    VaspBandsWorkchain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.bands.VaspBandsWorkchain>`"
    VaspRelaxWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.relax.VaspRelaxWorkChain>`"
    VaspConvergenceWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.converge.VaspConvergenceWorkChain>`"
    calcfunction: "{py:class}`calcfunction <aiida.engine.calcfunction>`"
    workfunction: "{py:class}`calcfunction <aiida.engine.workfunction>`"
---

(using-workflows)=

# Workflows

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

```{toctree}
./inputs
./design_principles
./bundled
./writing_workflows
```


[vasp]: https://www.vasp.at
