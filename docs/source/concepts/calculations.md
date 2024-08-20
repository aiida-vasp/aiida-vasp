(calculations)=

# Calculations

We call an [AiiDA] triggered execution of [VASP] a [Calculation]. The calculation has to take some given input on [AiiDA] form and projects it to [VASP] in an understandable manner. Similarly, on termination, the calculation has to parse the [VASP] output files such that the results are understandable for [AiiDA] and possibly also user friendly.

A [Calculation] is a special class in [AiiDA] that is derived from a [Process] class which handles these (and more subtle) tasks.

There are three types of [Calculation] classes for [VASP], an all-purpose one for pure [VASP] calculations ({ref}`vasp_calculation`) and two specialized ones for using the [VASP interface to Wannier90] ({ref}`wannier_calculation`) and immigrating [VASP] calculations that have been performed outside AiiDA-VASP. A AiiDA-VASP calculation can be accessed by loading it using the {py:class}`CalculationFactory<aiida.plugins.factories.CalculationFactory>` from [AiiDA]:

```
$ some_calc = CalculationFactory('<plugin_namespace>.<calculation_name>')
```

from the `verdi shell`. If you want to load it from a python script, please have a look at [verdi_shell]. The `<plugin_namespace>` is always `vasp` for the AiiDA-VASP plugin. The `<calculation_name>` is the name of the file containing the module. For instance, for the {ref}`vasp_calculation` we would issue:

```
$ vasp_calc = CalculationFactory('vasp.vasp')
```

Calculations should be placed in the `aiida_vasp/calcs` folder.

The general user should not care too much about the [Calculation] itself as we believe it is better for the user to interact with [VASP], or the other calculators from the {ref}`workchains`.

[aiida]: https://www.aiida.net
[calculation]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/calculations.html
[process]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/processes.html
[vasp]: https://www.vasp.at
[vasp interface to wannier90]: https://cms.mpi.univie.ac.at/wiki/index.php/LWANNIER90
[verdi_shell]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/working_with_aiida/scripting.html
