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
    VaspInputGenerator: "{py:class}`VaspInputGenerator <aiida_vasp.protocols.generator.VaspInputGenerator>`"
    PresetConfig: "{py:class}`PresetConfig <aiida_vasp.protocols.generator.PresetConfig>`"
---

(workflow_inputs)=

# Setting inputs of a workflow

This section will provide a brief overview of the internals of the VASP workflows.

The input and outputs of the workflows as implemented as `WorkChain` are AiiDA's `Data` types.
The `Data` is a subclass of the `Node` class, which represents data that is stored in the database
and could be used by other  `Process` nodes.
A `WorkChain` has a set of pre-defined input and output ports (which can be dynamic, if needed) that
specifies the types of data that can be passed to and from it.

Some python native types (`float`, `dict`, `str`) have their `Data` counterparts, such as `Float`, `Dict`, `Str` - they can be used as inputs to the workflows directly, but the conversion still takes
place internally.

There are two ways to pass inputs to the workflows. The most general way is to pass a `dict` object
contains key-values pairs of the data to be passed to each input port of the workchain.

```python
from aiida.engine import run_get_node
from aiida.plugins import WorkflowFactory
workflow = WorkflowFactory('core.arithmetic.multiply_add')
node = run_get_node(workflow, x=Int(1), y=Int(2), z=Int(3)).node
print(node.outputs.result)  # 9
```

:::{note}
The first argument should be the workchain class, followed by keyword inputs for each input port.
The `run_get_node` function launches the workchain with the current python interpreter, and in
production environments one typically uses the `submit` function instead.
In this case the workchain is stored in the database and marked to be executed by the daemon.
:::

For more complex workflows, we typically construct a dictionary and use the `**inputs` syntax to pass it to function that launches the workchain.

## The `ProcessBuilder` class

The approach above is very general but can be cumbersome for complex workflows with many inputs.
In addition, the user must somehow *remember* all the input port names and their types.
To address this problem, AiiDA provides the `ProcessBuilder` class, which can be used to construct
the inputs for a workflow in a more structured and interactive way.
For example (to be run inside a `verdi shell`)

```python
from aiida.engine import run_get_node
from aiida.plugins import WorkflowFactory

builder = WorkflowFactory('vasp.v2.vasp').get_builder()
builder.parameters = Dict(dict={'incar': {'encut': 500, 'ismear': 0}})
builder.kpoints_spacing = 0.05
```

The `builder` object has attributes corresponding to the input ports of the `VaspWorkChain`.
The conversion and validation of the inputs is done automatically when it is assigned to the attribute.

## The `InputGenerator` class

The `ProcessBuilder` class is a convenient way to construct inputs for a workflow, but one still
has to write inputs explicitly. To make it easier to construct inputs, the plugin provides the {{ VaspInputGenerator }} class.
As the name suggests, it is used to update the inputs of an underlying `ProcessBuilder` object.
The main advantage is that it allows the user to start from a predefined set of input values which
can be modified or added to.

There two kinds of pre-defined defaults that a {{ VaspInputGenerator }} can uses.
The first one is based on [protocols](../concepts/protocols.md#using-protocol-and-inputgenerator),
which defines base parameters for calculations and workflows.
The protocols may contain the default INCAR tags, the k-points spacing to be used and the
pseudopotential configurations as well as higher level control parameters for workchains.

The default `balanced` protocol is stored in the `<root>/src/protocols/` folder with the
following content:

:::{literalinclude} ../../../src/aiida_vasp/protocols/vasp.yaml
:::


While protocols default how parameters of each calculation and workchain are defined.
The  {{ PresetConfig }} offers control at a higher level - it records the default input set to be used as well as any code-specific overrides needed.

The default configuration is stored in the `<root>/src/aiida_vasp/protocols/presets/default.yaml` with the following content:

:::{literalinclude} ../../../src/aiida_vasp/protocols/presets/default.yaml
:::

This default preset file is used for tests and documentation examples with `mock-vasp@localhost` code.

Using the {{ VaspInputGenerator }} class is intended to simply the input construction process.
For example, to construct a `VaspWorkChain` with the default INCAR tags, k-points spacing and pseudopotential for a silicon structure (`si_node`), can be a simple as:

```python
from aiida_vasp.protocols.generator import VaspInputGenerator

upd = VaspInputGenerator()
upd.get_builder(structure=si_node, code='<my_code>@<computer>')
upd.submit()
```

When not using {{ VaspInputGenerator }}, the `get_builder_from_protocol` method of the workchain can be used to obtain the `ProcessBuilder` directly.
Alternatively,  the workchain will have to specified either through a multi-line mini script using the `ProcessBuilder` obtained with `get_builder` method or a large nested dictionary for complex workflows.

Nevertheless, one should still inspect the actual input passed to the workchain, this can be done
by simply returning the `builder` attribute of the {{ VaspInputGenerator }} object.


```python
upd.builder  # Should print the input to each port namespace of the workchain
```

There are `InputGenerator` class specific to each class:

|   WorkChain class | InputGenerator class |
| ------------------- | ---------------------- |
| {py:class}`aiida_vasp.workchains.v2.relax:VaspRelaxWorkChain` | {py:class}`aiida_vasp.protocols.generator:VaspRelaxInputGenerator` |
| {py:class}`aiida_vasp.workchains.v2.vasp:VaspWorkChain` | {py:class}`aiida_vasp.protocols.generator:VaspInputGenerator` |
| {py:class}`aiida_vasp.workchains.v2.band:VaspBandsWorkChain` | {py:class}`aiida_vasp.protocols.generator:VaspBandsInputGenerator` |
| {py:class}`aiida_vasp.workchains.v2.band:VaspHybridBandsWorkChain` | {py:class}`aiida_vasp.protocols.generator:VaspHybridBandsInputGenerator` |
| {py:class}`aiida_vasp.workchains.v2.converge:VaspConvergenceWorkChain` | {py:class}`aiida_vasp.protocols.generator:VaspConvegenceInputGenerator` |


## Customization

In practice, one may want to  have their own default inputs
This can be achieved by creating a new `<preset_name>.yaml` file inside `~/.aiida-vasp/presets` with the desired settings. The default configuration shown above can be used as a starting point.

It is also possible to have your own **protocol** - simply place the YAML files in the same `~/.aiida-vasp/<workchain tag>/<alias>.yaml` folder.

The `<workchain tag>` is an alias for the workchains.
| WorkChain Tag | Workchain  Class |
| -------------- | ----------- |
| vasp | {py:class}`aiida_vasp.workchains.v2.vasp:VaspWorkChain` |
| relax | {py:class}`aiida_vasp.workchains.v2.relax:VaspRelaxWorkChain` |
| band | {py:class}`aiida_vasp.workchains.v2.relax:VaspBandsWorkChain` |
| band | {py:class}`aiida_vasp.workchains.v2.relax:VaspHybridBandsWorkChain` |
| conv | {py:class}`aiida_vasp.workchains.v2.converge:VaspConvergenceWorkChain` |

The `<alias>` is an user defined alias for the protocol set.

:::{note}
:::
One should be careful when modifying or extending existing **preset** or **protocol** files as they may render calculations results incompatible for comparison.
Thankfully, AiiDA still preservers the full provenance of the calculation can be traced as the actual inputs are faithfully stored in the database.
