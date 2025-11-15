# How to analyze results

## How to check the status of a calculation/workflow?

* Using `verdi process status <process_pk>`: this will display a tree-like diagram containing the called processes.
* For a `VaspCalculation`, one can use `verdi calcjob <sub_command> <calculation_pk>` commands to show its info, commonly used sub-commands are:
   *  `inputls`: List the input files.
   *  `inputcat`: Print an input file. The name of the file needs to be passed following the pk, if no default calculation input file is defined (default is `INCAR`). The submission script can be displayed by passing `_aiidasubmit.sh`.
   *  `outputls`: List the output files.
   *  `outputcat`: Same as `inputcat` but print an output file instead. The default output is `OUTCAR`.
   *  `remotecat`: Same a `outputcat` but can be used for running calculations.
   *  `gotocomputer`: This command will take the current shell to the running folder of submitted calculation, which is very useful for inspecting running calculation/check the correctness of the input files and the submission script.
   *  `res`: Print the results of a calculation to the screen. This will display the `misc` output of a `VaspCalculation`.
   * `cleanworkdir`: Clean the working directly of a calculation.

:::{note}
The `verdi calcjob` command is inspecting a `Calculation`. A workflow may launch many calculation. In this case, one can use `verdi process status` to find the *pk* of calculations that have been launched.
:::

* Finished workflows can be *dumped* to the disk using `verdi process dump` command. This will create a multi-level directory containing the launched processes.
* This plugin provides a command `aiida-vasp tools export` which can be used to export completed calculations and workchains. The output if similar to `verdi process dump` with some differences:
  * The input and output files of `VaspCalculation`s are collected into a single folder, mimicking normal VASP calculations.
  * A `--include-potcar` option can be passed so the `POTCAR` file of each calculation is re-created. This is not the case when using `verdi process dump`, since the exact `POTCAR` content is not included in the provenance graph in order for the data to be sharable (for licensing reasons).

## How to obtain forces and stress of each ionic step

By default, only the forces, stress and energies of the last ionic step are stored in the `misc` output.
If you want those for each ionic step, you can modify the parser to enable the output `trajectory` node:

```python
from aiida.orm import WorkflowFactory
VaspWorkChain = WorkflowFactory('vasp.vasp')
builder = VaspWorkChain.get_builder()
settings ={'parser_settings': {'include_node': ['trajectory']}}
builder.settings = settings
```
