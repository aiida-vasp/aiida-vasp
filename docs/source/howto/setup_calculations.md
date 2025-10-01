---
file_format: mystnb
kernelspec:
  display_name: Python 3
  name: python3

myst:
  substitutions:
    VaspWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.vasp.VaspWorkChain>`"
    VaspCalculation: "{py:class}`VaspCalculation <aiida_vasp.calcs.vasp.VaspCalculation>`"
    VaspBandsWorkchain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.bands.VaspBandsWorkchain>`"
    VaspRelaxWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.relax.VaspRelaxWorkChain>`"
    VaspConvergenceWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.converge.VaspConvergenceWorkChain>`"
    calcfunction: "{py:class}`calcfunction <aiida.engine.calcfunction>`"
    workfunction: "{py:class}`calcfunction <aiida.engine.workfunction>`"
    VaspInputGenerator: "{py:class}`VaspInputGenerator <aiida_vasp.protocols.generator.VaspInputGenerator>`"
    PresetConfig: "{py:class}`PresetConfig <aiida_vasp.protocols.generator.PresetConfig>`"
---

(workflow_inputs)=

# How to setup calculations

## How to know what inputs are allowed for a calculation/workchain

AiiDA workchains define their inputs using *input ports* and *port namespaces*. Each workchain exposes a set of input ports (such as `structure`, `parameters`, `kpoints`, etc.), and these can be grouped into namespaces for logical organization (e.g., `parameters.incar`). This structure allows for flexible and hierarchical input definitions, making it easier to manage complex workflows.

The easiest way to explore what inputs a workchain can take is to use the `ProcessBuilder` with tab completion

```python
from aiida.plugin import WorkFlowFactory
builder = WorkflowFactory('vasp.vasp').get_builder()
builder.<tab>
```

Alternatively, one can use `verdi` commandline tool to inspect a workchain:

```bash
verdi plugin list aiida.workflows vasp.vasp
```

The third way is to look into the source code, for example, the  `VaspWorkChain` has the following lines :

```python
class VaspWorkChain:

    @classmethod
    def define(cls, spec: ProcessSpec) -> None:
        super(VaspWorkChain, cls).define(spec)
        spec.expose_inputs(cls._process_class, exclude=('metadata',))
        spec.expose_inputs(
            cls._process_class, namespace='calc', include=('metadata',), namespace_options={'populate_defaults': True}
        )

        # Use a custom validator for backward compatibility
        # This needs to be removed in the next major release/formalized workchain interface
        spec.inputs.validator = validate_calc_job_custom
        spec.inputs['calc']['metadata']['options']['resources']._required = False

        spec.input('kpoints', valid_type=orm.KpointsData, required=False)
        spec.input(
            'potential_family',
            valid_type=orm.Str,
            required=True,
            serializer=to_aiida_type,
            validator=potential_family_validator,
        )
```

In the code above, the `spec.input` call define a series of ports that a calculation may take as inputs.
An input port may contain certain default value, and it may or may not be *required* by a calculation.
A more advanced usage is the `spec.expose_inputs` call, which **expose** existing input ports of another calculation/workchain to the current workchain.
In above, the inputs of a {{ VaspCalculation }} is exposed at the top level as well as nested in a `calc`.
However, the latter only contains a `metadata` port, which is a special input port that allow defining *options* with request resource and wall-time limits or a  {{ VaspCalculation }}.

:::{admonition} Recommended workflow design pattern
:class: dropdown

This design pattern is due to historical reasons. For new projects, we recommend exposing inputs of a sub-workchain/calculation in full inside a nested namespace instead by default.

For example, the following code exposes *all* ports of a {{ VaspCalculation }} except the `kpoints` port:

```python
class UserWorkChain:

    @classmethod
    def define(cls, spec: ProcessSpec) -> None:
        super(VaspWorkChain, cls).define(spec)
        spec.expose_inputs(
            cls._process_class, namespace='calc',
            exclude=('kpoints',)
        )

        spec.input('kpoints', valid_type=orm.KpointsData, required=False)
        spec.input('kpoints_spacing', valid_type=orm.Float, required=False)
        spec.input(
            'potential_family',
            valid_type=orm.Str,
            required=True,
            serializer=to_aiida_type,
            validator=potential_family_validator,
        )
```

The `kpoints` port is exposed at the top level, but when `self.exposed_inputs` is called with `agglomerate=True` (default), the parent namspace is searched with the defined ports.
Hence the  `kpoints` port at the top-level
will be gathered as the input for a {{ VaspCalculation }} automatically.

As an example, if one defines the inputs to a `UserWorkChain` like:
```
UserWorkChain
|- kpoints
|- potential_family
|- calc
    |-parameters
```

Calling `UserWorkChain.exposed_inputs(VaspWorkChain, 'calc', agglomerate=True)` will have the following
inputs gathered and ready to be passed to a {{ VaspCalculation }}

```
VaspCalculation
|- kpoints
|- potentials
|- parameters
|- ...
```

:::

## How to setup inputs of a calculation or a workflow

In this document we will learn how to pass necessary input to the calculation and workflows provided by aiida-vasp.

The input and outputs of the workflows as implemented as `WorkChain` are AiiDA's `Data` types.
The `Data` is a subclass of the `Node` class, which represents data that is stored in the database
and could be used by other  `Process` nodes.
A `WorkChain` has a set of pre-defined input and output ports (which can be dynamic, if needed) that
specifies the types of data that can be passed to and from it.

Some python native types (`float`, `dict`, `str`) have their `Data` counterparts, such as `Float`, `Dict`, `Str` - they can be used as inputs to the workflows directly, but the conversion still takes
place internally.

There are tree ways to pass inputs to the workflows. The most general way is to pass a `dict` object
contains key-values pairs of the data to be passed to each input port of the workchain.

### Using a `dict` object as inputs to a `Process`

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

:::{note}
The first argument should be the workchain class, followed by keyword inputs for each input port.
The `run_get_node` function launches the workchain with the current python interpreter, and in
production environments one typically uses the `submit` function instead.
In this case the workchain is stored in the database and marked to be executed by the daemon.
:::

For more complex workflows, we typically construct a dictionary and use the `**inputs` syntax to pass it to function that launches the workchain.

### The `ProcessBuilder` class

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

### The `InputGenerator` class

While `ProcessBuilder` class is a convenient, one still
has to write inputs explicitly.
To make it even easier to construct inputs, we provide the {{ VaspInputGenerator }} class,
which can be used to generate a `ProcessBuilder` object using pre-defined [protocols](../concepts/protocols)
The main advantage is that it allows the user to start from a predefined set of input values which
can be modified or added to.

There two kinds of pre-defined defaults that a {{ VaspInputGenerator }} can uses.
The first one is based on [protocols](../concepts/protocols),
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


### Customize protocols and presets

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

For example, to have a custom relaxation protocol, create a file at `~/.aiida-vasp/protocols/relax/custom.yaml` with the follow content:

```yaml
# Default input values for the workflow
default_inputs:
  verbose: False               # Verbosity of the workflow output
  base_workchain_protocol: balanced

# Protocol definitions
protocols:
  posonly:
    description: |
      A protocol for relaxation that only relax positions
    relax_settings:
      shape: false
      volume: false
      positions: true
  fixxy:
    description: "Relax only the z axis"
    vasp:
      parameters:
        ioptcell: 0 0 0 0 0 0 0 0 1
    base_workchain_protocol: balanced
```

The protocol above can be referenced using `VaspRelaxWorkChain.get_builder_from_protocol(...., protocol="posonly@custom"`


:::{caution}
One should be careful when modifying or extending existing **preset** or **protocol** files as they may render calculations results incompatible for comparison.
Thankfully, AiiDA still preservers the full provenance of the calculation can be traced as the actual inputs are faithfully stored in the database.
:::

## How to fix the atoms during relaxation

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

:::{caution}
The `T` and `F` applies to the **direct** (fractional) coordinates.
To fix the **cartesian** coordinates, the $$lattice vectors needs to align with the x, y, z direction respectively.
:::


## How set initial magnetization for magnetic calculations

VASP uses the `MAGMOM` tag in the INCAR file for the initial magnetic momenets.
This tag can be an explicit list of moments for each specie:

```
MAGMOM = 0.0 1.0 5.0 5.0
```

or a compact format with `*` indicating repeats:

```
MAGMOM = 0.1*1 1.0*1 5.0*2
```

Very often, we want the initial magnetic moments to be specie-dependent.

We can of course explictly set the `MAGMOM` tags programatically::

```
config = {'O': 0.0, 'Fe': 5.0}
builder = VaspWorkChain.get_builder()
builder.structure = structure
builder.parameters = {
  'incar': {
    ....,
    'magmom': [config.get(site.kind_name, 0.6) for site in structure.sites]
  }
}
```

This will assign the initial magnetic moment of all `Fe` atoms to be 5.0 and `O` atoms to be 0.0, other species will have 0.6 as the initial magnetic moment.

In fact, {{ VaspWorkChain }} provides an input port to do exactly this.
The above code is equivalent to:

```python
config = {'O': 0.0, 'Fe': 5.0, 'default': 0.6}
builder = VaspWorkChain.get_builder()
builder.magmom_mapping = config
```

In addition, the `ISPIN` tag will be set to `2` if it is not explicitly defined in the input.

:::{note}
At the time of writing, `magmom_mapping` port cannot be used for setting the initial three-dimensional spin for none-collinear spin calculations.
:::

## How configure LDA+U calculations

LDA+U calcualtion in VAPS requires the following tags to be set in INCAR:

* `LDAU`: should be set to `.TRUE.` as an overall switch of +U
* `LDAUTYPE`: determines the type of the LDA+U formulism. Usually `2` is used which uses only `LDAUU`.
* `LDAUU`: sets the $$U$$ value of each specie
* `LDAUJ`: sets the $$J$$ value of each specie. The value of $$U-J$$ is often referred as $$U_{eff}$$ for a sepcific specie. This parameter is only used when `LDAUTYPE` is `1`.
* `LDAUL`: controls which angular momentum channel the +U should be applied for each specie. For *3d* transition metals, the angular momentum channel is `2`. The U is not applied if it is set to `-1`.

For example, to add U on Ni atoms with $$U_eff=6 \mathrm{eV}$$ in NiO, we can set

```
LDAU = .TRUE.
LDAUTYPE = 2
LDAUU = 6.0 0.0
LDAUL = 2 -1
```

This assumes that Ni comes first in the POSCAR and POTCAR, followed by O.

We can explicitly define these tags in the inputs to various workchain.
However, it is easier (and less prone to mistake) to use the `ldau_mapping` port for automatically generating these tags:

```python
config = {'mapping': {'Ni': [2, 6.0], 'Fe': [2, 4.0]}}
builder = VaspWorkChain.get_builder()
builder.structure = structure
builder.ldau_mapping = config
```

which sets $$U_{eff} = 6.0 \mathrm{eV}$$ for Ni's $$l=2$$ angular momentum channel, e.g. its $$3d$$
Internally, the {{ VaspWorkChain }} uses {py:func}`get_ldau_keys <aiida_vasp.utils.ldau:get_ldau_keys>` function to generate the INCAR tags and the available keys can be found in its docstring.
