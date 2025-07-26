---
file_format: mystnb
kernelspec:
  name: python3
myst:
  substitutions:
    VaspWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.vasp.VaspWorkChain>`"
    VaspCalculation: "{py:class}`VaspCalculation<aiida_vasp.calcs.vasp.VaspCalculation>`"
---

(parsing)=

# Parsing

:::{note}
This notebook can be downloaded as **{nb-download}`parsing.ipynb`** and {download}`parsing.md`
:::

AiiDA-VASP provides flexible parsing of VASP output files to store data in the AiiDA database and repository.

The user interface for configuring the parsing settings takes place in the `settings['parser_settings']` dictionary entry. The default `parser_settings` is presently:


```{code-cell}
from aiida_vasp.parsers.vasp import ParserSettingsConfig
ParserSettingsConfig().dict()
```

% Use executable cell to display the default parser_settings
%```{literalinclude} ../../../aiida_vasp/parsers/vasp.py
%:starts-after: class ParserSettingsConfig
%:ends-before: class VaspParser
%```

These settings can be configured by settings the `parser_settings` inside the `settings` input node to the {{VaspCalculation}} or {{VaspWorkChain}} processes.

The parser is responsible for extracting information for a VASP calculation and store them as {py:class}`Data <aiida.orm.Data>` nodes in the database.

## Error checking

The parser will check for any error detected for the underlying VASP calculation.
This information is stored in the `notification` quantity stored in the `misc` output node, which contains a list of error/warnings detected in the STDOUT of the VASP calculation.
By default, a non-zero exit state is returned if any critical error is found.
The default list of critical errors is defined under the `critical_notification_errors` key inside the `parser_settings`.
Additional settings may also be supplied under `parser_settings` to modify the behaviors:

- `ignore_notification_errors`: a boolean value to control whether all notifications should be ignored, defaults to `False`.
- `critical_notification_errors`: a list with items for controlling which errors present in the stdout to be regarded as `critical` or not.
- `check_completeness`: whether to check the completeness of the retrieved files, defaults to `True`.
- `critical_objects`: a list of objects (files) that must be present as the outputs of a calculation.
- `check_errors`: whether to check for errors in the VASP output, defaults to `True`.

## Configuring the output nodes

The parser will generate a set of output nodes based on the parsed data.
Most of the common outputs are placed in the `misc` output node ({py:class}`Dict <aiida.orm.Dict>`),
and those that are typically large arrays are stored in the `array` output node ({py:class}`ArrayData <aiida.orm.ArrayData>`).

Other output nodes includes:

- `structure`: the final structure of the calculation ({py:class}`StructureData <aiida.orm.StructureData>`).
- `bands`: the electronic band structure ({py:class}`BandsData <aiida.orm.BandsData>`).
- `dos`: the density of states ({py:class}`ArrayData <aiida.orm.ArrayData>`).
- `born_charges`: the Born effective charges ({py:class}`ArrayData <aiida.orm.ArrayData>`).
- `dielectrics`: the dielectric function ({py:class}`ArrayData <aiida.orm.ArrayData>`).
- `hessian`: the Hessian matrix ({py:class}`ArrayData <aiida.orm.ArrayData>`).
- `dynmat`: the dynamical matrix ({py:class}`ArrayData <aiida.orm.ArrayData>`).
- `wavecar`: the wave function ({py:class}`ArrayData <aiida.orm.ArrayData>`).

Note that the parser will only generate output nodes for those data that are present in the retrieved files.
For simplicity, all available quantities are collected from the parser for each retrieved files.
One can turn on/off of these outputs by modifying the *included* quantities as well as turning on/off the node output directly.
For example, the `bands` node is essentially made of two quantities: `eigenvalues` and `occupancies`, but the
band structure output if off by default at the node level.
This is because most of the calculation do contain
these two quantities, but only those from a dedicated band structure calculation are meaningful in analysis.
We achieve this by including the `bands` node in the `DEFAULT_EXCLUDED_NODE`.
For a band structure calculations, one should set `settings['parser_settings']['include_node']` to `['bands']` to indicate that the `bands` node should be created.

Likewise, certain quantities/nodes are excluded by default:

```{code-cell}
from aiida_vasp.parsers.vasp import DEFAULT_EXCLUDED_NODE, DEFAULT_EXCLUDED_QUANTITIES
print(DEFAULT_EXCLUDED_NODE)
print(DEFAULT_EXCLUDED_QUANTITIES)
```
