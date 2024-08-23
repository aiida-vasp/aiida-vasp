---
substitutions:
  VaspWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.vasp.VaspWorkChain>`"
---
(workflows)=

# Workflows

By composing one or several [Workchain] classes, one can make a workflow. A workflow can of course only consist of one rather trivial operation, but is becoming truly useful when the task at hand is complicated and time consuming to construct, inspect, nourish and analyze by hand. With the introduction of high-throughput calculations and more composed single calculations using [VASP] such a feature is welcoming.

The bundled workflows in AiiDA-VASP is currently composed of geometry optimization (relax), convergence, band as well as the basic workchain  ({{ VaspWorkChain }}).

We would like to encourage users to build workchains and/or compose existing ones into more advanced workflows that we can all share and benefit from.

One should note that the advantage of using a provenance-preserving engine like AiiDA is that you do
not have to define a workflow in order to have the calculations steps recorded and stored.
It is perfectly fine to conduct exploration studies using the basic workchains and use {py:func}`calcfunction <aiida.engine.processes.functions.calcfunction>`  to link the outputs/inputs together for provenance.

[vasp]: https://www.vasp.at
[workchain]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/workflows.html#work-chains
