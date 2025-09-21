---
myst:
  substitutions:
      VaspWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.vasp.VaspWorkChain>`"
---

(parameters)=

# Protocols

Protocols in `aiida-vasp` are YAML files (see `src/aiida_vasp/protocols/`) that specify recommended or standardized sets of input parameters for different types of calculations (e.g., relaxation, band structure, convergence). Each protocol can define default values for any input port or namespace, and can also provide multiple named protocols (such as `fast`, `balanced`, `stringent`) for different accuracy or speed requirements.

:::{hint}
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
