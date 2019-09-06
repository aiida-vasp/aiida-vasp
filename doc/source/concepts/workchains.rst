.. _workchains:

==========
Workchains
==========

`AiiDA`_ utilizes the concept of a `Workchain`_ as a building block to creating more complex workflows. AiiDA-VASP honours this. The idea is that the workchains should be modular such that one can pick what is neaded and build a workflow that suits the problem at hand. Similar to `Calculation`_, the `Workchain`_ is a special class in `AiiDA`_ that is derived from a `Process`_ class.

A `Workchain`_ can be loaded by utilizing the :py:class:`aiida.plugins.WorkflowFactory`::

  $ some_workchain = WorkflowFactory(<plugin_namespace>.<workchain_name>)

from the ``verdi shell``. If you want to load it from a python script, please have a look at `verdi_shell`_. The ``<plugin_namespace>`` is always ``vasp`` for the AiiDA-VASP plugin. The ``<workchain_name>`` is the name of the file containing the module. For instance, for the :ref:`VASP workchain` we would issue::

  $ vasp_workchain = WorkflowFactory(vasp.vasp)

Workchains should be placed in the ``aiida_vasp/workchains`` folder.

AiiDA-VASP is delivered with several preconstructed workchains to perform dedicated tasks and we believe these, instead the the :ref:`calculations` should be the main entry point for users. They are currently:

The VASP workchain
------------------
Performs the low level interactions with the :ref:`vasp_calculation` and is accesible by utilizing the :py:class:`aiida.plugins.WorkflowFactory` with the key ``vasp.vasp``. For additional details on how to interact with this workchain, see :ref:`vasp_workchain`.

The verify workchain
--------------------
Constructed to handle input/output verifications that is not related to numerics or the construction of combinations of workchains. Calls :ref:`vasp_workchain`. Currently, this workchain is not doing anything. For additional details on how to interact with this workchain, see :ref:`verify_workchain`.

The relaxation workchain
------------------------
Handles the process of relaxing the structure and calls :ref:`verify_workchain`. For additional details on how to interact with this workchain, see :ref:`relax_workchain`.

The convergence workchain
-------------------------
Checks the k-point grid and plane wave cutoff for convergence on selected parameters. Calls :ref:`relax_workchain`. For additional details on how to interact with this workchain, see :ref:`converge_workchain`. We consider this the main entrypoint for a regular calculation with `VASP`_.

The bands workchain
-------------------
Enables the calculation of electronic band structures by standardizing the selected paths in reciprocal space using `SeeKpath`_. Calls the :ref:`vasp_workchain`. For additional details on how to interact with this workchain, see :ref:`bands_workchain`. If a user wants to calculate the electronic band structure they should use the :ref:`master_workchain` as the main entry point.

The master workchain
--------------------
The idea of this workchain is to ultimately be the main entry point, such that a user can select what properties they want calculated and the master workchain then composes a workflow to enable such extraction. Currently only the calculation of the electronic band structure is enabled. But this serves as a nice introductory example that can be easily expanded. Calls any relevant workchain, depending on the chosen input parameters. For additional details on how to interact with this workchain, see :ref:`master_workchain`.

.. _AiiDA: https://www.aiida.net
.. _Workchain: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/workflows.html#work-chains
.. _Process: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/processes.html
.. _Calculation: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/calculations.html
.. _VASP: https://www.vasp.at
.. _`SeeK-path`: https://github.com/giovannipizzi/seekpath
.. _verdi shell: https://aiida.readthedocs.io/projects/aiida-core/en/latest/working_with_aiida/scripting.html
