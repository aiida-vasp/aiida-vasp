---
myst:
  substitutions:
      VaspWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.vasp.VaspWorkChain>`"
---

(parameters)=

# Passing parameters to VASP

Before describing how parameter passing works in this plugin it is worthwhile to restate that the design principle is that all higher lying workchains ultimately call the {{ VaspWorkChain }}
which should handle [VASP] specific translations and setups in order to execute your problem with [VASP]. At that point what we in general call parameters are fully converted to INCAR tags or flags in POSCAR, for instance in the case of selective dynamics.

:::{note}
In this documentation, there is the parameters, which is the general description of something you can adjust to get some specific behavior, or `parameters` which is
a dedicated input parameter.
:::

We now describe how parameters can be passed in the plugin. We separate between passing parameters directly to the `VaspCalculation` ({ref}`vasp_calculation`), the {py:class}`VaspWorkChain<aiida_vasp.workchains.vasp.VaspWorkChain>` (or any workchain ultimately calling {py:class}`VaspWorkChain<aiida_vasp.workchains.vasp.VaspWorkChain>`).
The latter being the recommended approach, unless you have very specific use-cases that warrants interacting with the {py:class}`VaspCalculation<aiida_vasp.calcs.vasp.VaspCalculation>`.

## Direct to VASP calculations

This is the least used approach. Defining inputs of `VaspCalculation` requires explicitly setting all
relevant inputs just like defining the calculations via input fields (manual "crafted" calculations).
This is by design as we want to fully capture the provenance of each calculation and ensure reproducibility.

The INCAR tags are directly defined under the `parameters` input node as a `orm.Dict` object.
These tags should be in lower case by convention.

# Using `protocol` and `InputGenerator`

## How inputs to a workchain are defined

AiiDA workchains define their inputs using *input ports* and *port namespaces*. Each workchain exposes a set of input ports (such as `structure`, `parameters`, `kpoints`, etc.), and these can be grouped into namespaces for logical organization (e.g., `parameters.incar`). This structure allows for flexible and hierarchical input definitions, making it easier to manage complex workflows.

The easiest way to explore what inputs a workchain can take is to use the `ProcessBuilder` with tab completion

```python
from aiida.plugin import WorkFlowFactory
builder = WorkflowFactory('vasp.vasp').get_builder()
builder.<tab>
```

Alternatively, one can use `verdi` commandline tool to inspect a workchain:

```base
verdi plugin list aiida.workflows vasp.vasp
```

The third way is to look into the source code, for example, the  `VaspWorkChain` has the following lines :

```python
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
In above, the inputs of a `VaspCalculation` is exposed at the top level as well as nested in a `calc`.
However, the latter only contains a `metadata` port, which is a special input port that allow defining *options* with request resource and wall-time limits or a `VaspCalculation`

:::{info}
This design pattern is due to historical reasons. For new projects, we recommend exposing inputs of a sub-workchain/calculation in full inside a nested namespace instead
:::

## Design of protocols

Protocols in `aiida-vasp` are YAML files (see `src/aiida_vasp/protocols/`) that specify recommended or standardized sets of input parameters for different types of calculations (e.g., relaxation, band structure, convergence). Each protocol can define default values for any input port or namespace, and can also provide multiple named protocols (such as `fast`, `balanced`, `stringent`) for different accuracy or speed requirements.

:::{info}

The default protocol `balanced` is equivalent to the `UCLRelaxSet` which uses a PBEsol functional. A slight difference is that LDA+U is not automatically applied while `UCLRelaxSet` applies U=4.0 eV for Ti and Fe by default.
:::

Protocols can be extended or overridden by users, either by providing their own YAML files (e.g., in `~/.aiida-vasp/protocols/`) or by passing overrides at runtime. The protocol system supports merging of nested dictionaries, so only the relevant parts need to be overridden.

The `ProtocolMixin` class provides methods to list available protocols, load protocol files, and generate the full set of inputs for a given protocol, including applying user overrides.


Below is an example of a protocol YAML file (e.g., for a relaxation workflow), with comments explaining each key:

```yaml
# Default input values for the workflow
default_inputs:
  relax_settings: # Settings for the relaxation (see RelaxOptions for all valid keys)
    shape: true                # Allow cell shape to change during relaxation
    volume: true               # Allow cell volume to change
    positions: true            # Allow atomic positions to relax
    steps: 60                  # Maximum number of ionic steps
    algo: cg                   # Relaxation algorithm (e.g., 'cg' for conjugate gradient)
    force_cutoff: 0.03         # Convergence threshold for forces (in eV/Å)
    convergence_on: true       # Enable convergence checks
    convergence_max_iterations: 5  # Max number of convergence cycles
    convergence_volume: 0.01   # Convergence threshold for volume change (fractional)
  verbose: False               # Verbosity of the workflow output

# The default protocol to use if none is specified
default_protocol: balanced

# Protocol definitions
protocols:
  balanced:
    description: |
      A balanced protocol for relaxation. See vasp.yaml for exact settings used.
  fast:
    relax_settings:
      convergence_volume: 0.05   # Looser volume convergence for faster runs
      force_cutoff: 0.05         # Looser force convergence for faster runs
```

The  `default_inputs` key defines the default input to the workflow.
In example, the content `relax_settings` port is explicitly defined, and will be converted to a `Dict` node automatically.
The `default_protcol` key defines the default protocol to use.
The `protocol` field defines a series of protocols as modifications to the default input.
The inputs to the workchain is constructed by merging the default input with the modifications.

It is worth noting that the VASP incar tags are missing in the example above.
This is because by default the protocols are cascaded to the lower level workchains.
In other words, the `balanced` protocol of `VaspWorkChain ` is applied automatically.


## Writing your own protocol


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

## Using `InputGenerator` classes to launch workchains

To simplify launching workchains with protocol-based inputs, `aiida-vasp` provides a set of `InputGenerator` classes (see `src/aiida_vasp/protocols/generator.py`). These classes automate the process of:

- Loading protocol and preset files
- Generating a builder with all required inputs set according to the chosen protocol
- Allowing interactive updates to common parameters (e.g., INCAR tags, computational resources, options)
- Supporting code-specific and user-specific presets

**Example usage:**

```python
from aiida_vasp.protocols.generator import VaspRelaxInputGenerator

# Instantiate the generator with a preset and protocol
gen = VaspRelaxInputGenerator(preset_name='default', protocol='balanced')

# Generate a builder for a given structure and code
builder = gen.get_builder(structure=my_structure, code='vasp@mycluster')

# Optionally update INCAR tags or computational options
gen.set_incar({'encut': 520, 'ismear': 0})
gen.set_options({'max_wallclock_seconds': 3600})

# Submit the workchain
from aiida.engine import submit
submit(builder)
```

This approach ensures that all inputs are consistent with the selected protocol, while still allowing for easy customization. Protocols and presets can be listed and inspected programmatically, and users can create their own YAML files to define custom protocols for specific projects or systems.



# Using `Inputset` and `BuilderUpdater` (deprecated)

These are old approaches originally designed for simple personal projects, but later merged into `aiida-vasp`.
The documentation shown here may be out-dated

## Using `VaspWorkChain`

At the first glance, the `VaspWorkChain` is just like a `VaspCalculations` but they are different in
several aspects. We do not go into the details here.
However, for a `VaspWorkChain` the `parameters` input may contain human-readable key-value pairs defining
how the INCAR tags should be set.
To set the INCAR tags directly, simply define the key-value pairs in the `incar` namespace of the `parameters` input node.
The workchain will workout the actual INCAR tags to be used and pass them to `VaspCalculation`.
In addition, the user may supply `potential_family` and `potential_mapping` to a `VaspWorkChain` for
defining the POTCAR files to be used.
There are a few other inputs such as `ldau_mapping`, `kpoints_spacing` that can be set.

## Using `VaspBuilderUpdater`

This is the easiest and recommended way to construct workflows as the inputs are automatically
constructed from presets that are stored as files.
The user may define their own custom inputs preset for specific projects, and the only input
required is the `structure`.

For example:

```python
from aiida_vasp.common import VaspBuilderUpdater

>>> upd = VaspBuilderUpdater("MyInputPreset").apply_preset(structure, label='My Awesome Calculation')
>>> upd.builder   # Inspect the builder - alway good to check if everything is as expected
>>> upd.submit()   # Submit the calculation to the daemon
```

## Other workchains

Some workchains may have their own specific parameters, for example, the `relax_settings` input for
a `VaspRelaxWorkChain` or the `band_settings` input for a `VaspBandWorkChain`. These parameters are
controls how the workchain behaves.

The convention of these workchains is to have the `structure` input and other settings
in the root namespace, and the other inputs (typically that of the `VaspWorkChain` inside the `vasp` namespace).
This way, higher level workchain can be defined easier by just exposing the relavant inputs of the
lower level workchains.

## Allowing custom [VASP] tags

In case you for instance perform developments in the [VASP] code, sometimes it makes sense to add a new [VASP] tag. This can be supplied in `settings.unsupported_parameters` as dict with the following specifications:

```
unsupported_parameters = {'my_unsupported_parameters': {
'default': 1.0,
'description': 'Some description',
'type': float,
'values': [1.0, 2.0]
}}
builder.settings = Dict(dict={'unsupported_parameters': unsupported_parameters})
```

Alternatively, the validation can be turned off entirely by setting `skip_parameters_validation` to `True` under `settings`, for example:

```
builder.settings = Dict(dict={'skip_parameters_validation': True})
```

The above works for both {py:class}`VaspWorkChain<aiida_vasp.workchains.vasp.VaspWorkChain>` and {py:class}`VaspCalculation<aiida_vasp.calcs.vasp.VaspCalculation>`.
In the latter case, if any of `skip_parameters_validation` or `unsupported_parameters` are present in the `settings` input node, the validation is turned off completely.

[vasp]: https://www.vasp.at
