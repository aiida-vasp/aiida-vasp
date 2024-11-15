---
myst:
  substitutions:
    VaspWorkChain: "{py:class}`VaspWorkChain <aiida_vasp.workchains.v2.vasp.VaspWorkChain>`"
    calcfunction: "{py:class}`calcfunction <aiida.engine.calcfunction>`"
    workfunction: "{py:class}`calcfunction <aiida.engine.workfunction>`"
---
(workflows)=

# Workflows

The [Workchain] class is the central piece that enables workflows to be run with aiida-vasp.
By composing one or several [Workchain] classes, one can make a workflow.
As single [WorkChain] class may launch one or several calculations, or it may launch children [WorkChain]s to achieve the designed functionality.

For  any short-running python code, the workchain can run them directly as {{ calcfunction }} or {{ workfunction }} directly, and the provenance will be recorded accordingly.

It is important to note that however, long-running computational *should not* but run directly in the code that it
will delay or block the operation of the [daemon].

We would like to encourage users to build workchains and/or compose existing ones into more advanced workflows that we can all share and benefit from.
You may want to visit [this page](https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/write_workflows.html) to learn more about WorkChains and how to build them.

One should note that the advantage of using a provenance-preserving engine like AiiDA is that you do
not have to define a workflow in order to have the calculations steps recorded and stored.
It is perfectly fine to conduct exploration studies using the basic workchains and use {py:func}`calcfunction <aiida.engine.processes.functions.calcfunction>`  to link the outputs/inputs together for provenance.

[vasp]: https://www.vasp.at
[workchain]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/workflows.html#work-chains
[daemon]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/topics/daemon.html
