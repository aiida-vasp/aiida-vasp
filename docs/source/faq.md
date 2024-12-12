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
