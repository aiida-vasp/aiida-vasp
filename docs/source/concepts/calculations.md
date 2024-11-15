---
myst:
  substitutions:
    VaspNEBCalculation: "{py:class}`VaspNEBCalculation <aiida_vasp.calcs.neb.VaspNEBCalculation>`"
    VaspCalculation: "{py:class}`VaspCalculation<aiida_vasp.calcs.vasp.VaspCalculation>`"
---

(vasp_calculation)=

# Calculations

We call an AiiDA triggered execution of VASP a [calculation job].
The calculation job has to take some given input on AiiDA form and projects it to VASP in an understandable manner.
Similarly, on termination, the calculation job has to parse the VASP output files such that the results are understandable for AiiDA and possibly also user friendly.

The calculation job is represented by a `CalcJobNode` in the database, which is a special class in AiiDA that is derived from a ProcessNode class which handles these (and more subtle) tasks.

There are two types of calculation classes for VASP, an all-purpose one for pure VASP calculations {{ VaspCalculation}} and the other is for NEB calculations with VTST tools {{ VaspNEBCalculation }}.

A AiiDA-VASP calculation can be accessed by loading it using the {py:class}`CalculationFactory<aiida.plugins.factories.CalculationFactory>` from [AiiDA]:

```
$ some_calc = CalculationFactory('<plugin_namespace>.<calculation_name>')
```

from the `verdi shell`. If you want to load it from a python script, please have a look at [verdi_shell]. The `<plugin_namespace>` is always `vasp` for the AiiDA-VASP plugin. The `<calculation_name>` is the name of the file containing the module. For instance, to load the two calculations mentioned above,  we would issue:

```
$ vasp_calc = CalculationFactory('vasp.vasp')
$ vasp_neb = CalculationFactory('vasp.neb')
```

Calculations should be placed in the `src/aiida_vasp/calcs` folder.

The general user should not care too much about the calculation itself as we believe it is better for the user to interact with VASP, or the other calculators from the workchains.
Nevertheless, defining calculations is a crucial part of the plugin as it defines how we interact with VASP and how the input and output data are stored in the database.

[verdi_shell]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/interact.html#how-to-interact-scripts
