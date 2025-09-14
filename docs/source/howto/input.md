# Defining calculation inputs

## Ways to define calculation inputs


The calculation and workchains are *processes* in the terminology of AiiDA.
Each are configured via a series of input **nodes** which are {py:class}`Data <aiida.orm:Data>` object.

The two common ways to launch a *process* is by using the `submit` or `run_get_node` functions in `aiida.engine`.

```python
from aiida.engine import submit

submit(Process, **inputs)
```

where `Process` is a class of the process to be launched. In aiida-vasp, it may be `VaspCalculation`, `VaspWorkChain` or other provided processes.
The `inputs` is a dictionary containing a nested key-value pairs defining inputs for each port of the process.
A typical `inputs` dictionary for `VaspWorkChain` looks like

```python
inputs = {
  'structure': si_structure,  # An instance of aiida.orm.StructureData
  'parameters': incar_tags,   # An instance of aiida.orm.Dict
  'calc':
    {'options':
      {'resources':
         {
          'num_machines': 1
         }
      }
    },
  # ....
}
```




## Where to pass X input for Y calculation/workflow?

The inputs to a workchain/calculation must be set to the correct input *ports*.
The easiest way is to use obtain a `builder` and use tab completion to browse through the ports:

```python
VaspWorkChain = WorkflowFactory("vasp.vasp")
builder = VaspWorkChain.get_builder()
builder.<tab>
```

Each port also has its owns documentation. If using IPython, you can use `?builder.<port>` to see its documentation.
If using Jupyter Notebook, `<ctrl> + <tab>` with the cursor at the end of the `<port>` should be able to trigger the documentation pop-out window.

Alternative, one can use `verdi plugin list aiida.workflows <entrypoint>` to print documentation of a workflow.
For example, to see the documentation of `VaspWorkChain`:

```
verdi plugin list aiida.workflows vasp.vasp
```

However, this does not alway display *expose inputs* correctly, so the `builder` method is preferred.


## Fixing the atoms during relaxation

Atoms may be fined by setting the `dynamics` input port using a dictionary:

```python
dynamics = {
  'positions_dof': [
    [1, 1, 1],
    [0, 0, 0],
    [0, 0, 1],
  ]
}
```

This means to:

```
T T T
F F F
T T F
```

in the generated POSCAR file for the three atoms for selective dynamics.
The `T` means that the atom is allow to move in this degree of freedom, and `F` means that the atom is fixed in this direction.

For example, if one wants to completely fix all atoms between  $\mathrm{2 \AA}$ to  $\mathrm{4 \AA}$ in the z direction:

```python
z = builder.structure.get_ase().positions[:, 2]
to_fix  = (z < 4) & (z > 2)
dof = [[1, 1, 1] if not fix else [0, 0, 0] for fix in to_fix]

builder.dynamics = {
  'positions_dof': dof
}
```

:::{warning}
The `T` and `F` applies to the **direct** (fractional) coorindates.
To fix the **cartesian** coorindates, the $$lattice vectors needs to align with the x, y, z direction respectively.
:::
