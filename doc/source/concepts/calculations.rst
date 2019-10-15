.. _calculations:

Calculations
============
We call an `AiiDA`_ triggered execution of `VASP`_ a `Calculation`_. The calculation has to take some given input on `AiiDA`_ form and projects it to `VASP`_ in an understandable manner. Similarly, on termination, the calculation has to parse the `VASP`_ output files such that the results are understandable for `AiiDA`_ and possibly also user friendly.

A `Calculation`_ is a special class in `AiiDA`_ that is derived from a `Process`_ class which handles these (and more subtle) tasks.

There are three types of `Calculation`_ classes for `VASP`_, an all-purpose one for pure `VASP`_ calculations (:ref:`vasp_calculation`) and two specialized ones for using the `VASP interface to Wannier90`_ (:ref:`wannier_calculation`) and immigrating `VASP`_ calculations that have been performed outside AiiDA-VASP (:ref:`immigrator_calculation`). A AiiDA-VASP calculation can be accessed by loading it using the :py:class:`aiida.plugins.CalculationFactory` from `AiiDA`_::

  $ some_calc = CalculationFactory('<plugin_namespace>.<calculation_name>')

from the ``verdi shell``. If you want to load it from a python script, please have a look at `verdi_shell`_. The ``<plugin_namespace>`` is always ``vasp`` for the AiiDA-VASP plugin. The ``<calculation_name>`` is the name of the file containing the module. For instance, for the :ref:`vasp_calculation` we would issue::

  $ vasp_calc = CalculationFactory('vasp.vasp')

Calculations should be placed in the ``aiida_vasp/calcs`` folder.

The general user should not care too much about the `Calculation`_ itself as we believe it is better for the user to interact with `VASP`_, or the other calculators from the :ref:`workchains`.

.. _Process: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/processes.html
.. _Calculation: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/calculations.html
.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _VASP interface to Wannier90: https://cms.mpi.univie.ac.at/wiki/index.php/LWANNIER90
.. _verdi_shell: https://aiida.readthedocs.io/projects/aiida-core/en/latest/working_with_aiida/scripting.html
