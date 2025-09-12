# Frequently Asked Questions (FAQ)

## How to obtain forces and stress of each ionic step?

By default, only the forces, stress and energies of the last ionic step are stored in the `misc` output.
If you want those for each ionic step, you can modify the parser to enable the output `trajectory` node:

```python
from aiida.orm import WorkflowFactory
VaspWorkChain = WorkflowFactory('vasp.vasp')
builder = VaspWorkChain.get_builder()
settings ={'parser_settings': {'include_node': ['trajectory']}}
builder.settings = settings
```

## Why the parser system is so complex?

Unfortunately a single DFT calculation generates lots of data, only some of which are useful in most cases (i.e. energy, forces, stress) and it depends on the type of calculations run.

A further complication is that VASP generate multiple output files and some *quantities* are repeated in different files.

The role of the parser is to parse *quantities* from the files and organize them into different output nodes in a meaningful way. In addition, to avoid overfloing the database, some quantities / output nodes are excluded by default, and the user can choose to include them by setting `'include_quantity'` or `'include_node'` in the `'parser_settings'` dictionary inside the `settings` input node of the `VaspCalculation`/`VaspWorkChain`.

The current logic of the parser system works like this:

- Instantiate all content parser for each kind of output file, which essentially parse everything and store them as their own attributes
- Collect quantities into a nested dictionary. All quantities declared as available by the content parser are collected unless they are explicitly _excluded_.
    - There is a default list of excluded quantities, and the user can include them by setting `include_quantity: ['<quantitity']` in parser settings
    - The main reason for having a list of default excluded quantities is because some quantities are not needed in most cases, but the *node* containing them should still be created (e.g. the `misc` output).
- Try to compose all nodes except those are excluded by default (again, can be overridden by `include_node: ['<node>']`).
    - If a node cannot be composed due to lack of quantity, we simply skip it, as it is the responsibility of the `CalcJob` and the higher-level workflows to check for required output.
    - Again the reason for having excluded nodes is that some nodes are only needed in specific cases, but the underlying quantities are always available. For example, the `'eigenvalues'` are avaliable in every single calculation, but they are mostly only needed for constructing the bands structure.
- Finally, we collect the compose nodes and store them under the `outputs` attribute which is a dictionary.

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

## How to check the results of a calculation/workflow?

* Using `verdi process status <process_pk>`: this will display a tree-like diagram containing the called processes.
* For a `VaspCalculation`, one can use `verdi calcjob <sub_command> <calculation_pk>` commends to show its info, commonly used sub-commends are:
   *  `inputls`: List the input files.
   *  `inputcat`: Print an input file. The name of the file needs to be passed following the pk, if no default calculation input file is defined (default is `INCAR`). The submission script can be displayed by passing `_aiidasubmit.sh`.
   *  `outputls`: List the output files.
   *  `outputcat`: Same as `inputcat` but print an output file instead. The default output is `OUTCAR`.
   *  `remotecat`: Same a `outputcat` but can be used for running calculations.
   *  `gotocomputer`: This command will take the current shell to the running folder of submitted calculation, which is very useful for inspecting running calculation/check the correctness of the input files and the submission script.
   *  `res`: Print the results of a calculation to the screen. This will display the `misc` output of a `VaspCalculation`.
   * `cleanworkdir`: Clean the working directly of a calculation.

:::{info}
The `verdi calcjob` command is inspecting a `Calculation`. A workflow may launch many calculation. In this case, one can use `verdi process status` to find the *pk* of calculations that have been launched.
:::

* Finished workflows can be *dumped* to the disk using `verdi process dump` command. This will created a multi-level directly containing the launched processes.
* This plugin provides a command `aiida-vasp tools export` which can be used to export completed calculations and workchains. The output if similar to `verdi process dump` with some differences:
  * The input and output files of `VaspCalculation`s are collected into a single folder, mimicking normal VASP calculations.
  * A `--include-potcar` option can be passed so the `POTCAR` file of each calculation is re-created. This is not the case when using `verdi process dump`, since the exact `POTCAR` content is not included in the provenance graph in order for the data to be sharable (for licensing reasons).

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




## I found the documentation to be out-dated and missing certain contents!


The documentation may not be update-to-date as new features are added to the code.
Sometimes it is because the developers feel the feature is not *completed* with time constraints.
In other cases, it may just because they forgot....

Feedbacks and suggestions on the documentation would be extremely helpful. Please create an issue if you have any in mind!
