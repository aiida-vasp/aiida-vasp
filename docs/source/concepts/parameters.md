---
myst:
  substitutions:
      VaspWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.vasp.VaspWorkChain>`"
---

(parameters)=

# Parameters

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
from aiida_vasp.workchains.v2.common import VaspBuilderUpdater

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
