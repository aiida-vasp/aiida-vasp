.. _vasp_workchain:

==============
VASP workchain
==============

Required inputs
^^^^^^^^^^^^^^^

The VaspWorkChain requires a number of inputs, these comprise the minimum set of information to run a `VASP`_ calculation from `AiiDA`_.

* ``code``: an `AiiDA`_ :py:class:`aiida.orm.nodes.data.Code` instance, describes the VASP executable and holds a reference to the ``Computer`` instance on which it lives.
* ``structure``: an `AiiDA`_ :py:class:`aiida.orm.nodes.data.StructureData` or :py:class:`aiida.orm.nodes.data.CifData` instance, describes the structure on which VASP is to be run.
* ``kpoints``: an `AiiDA`_ :py:class:`aiida.orm.nodes.data.KpointsData` instance, describing the kpoints mesh or path.
* ``potential_family``: an `AiiDA`_ :py:class:`aiida.orm.nodes.data.Str` instance, the name given to a set of uploaded POTCAR files.
* ``potential_mapping``: an `AiiDA`_ :py:class:`aiida.orm.nodes.data.Dict` instance, containing an entry for at least every kind name in the ``structure`` input with the full name of the POTCAR from the ``potential_family``. Example: ``{'In1': 'In_d', 'In2': 'In_h'}``.
* ``parameters``: an `AiiDA`_ :py:class:`aiida.orm.nodes.data.Dict` instance. Please consult the documentation on how parameters are handled (:ref `parameters`) for details, particularly the section pertaining to the ``VaspWorkChain``.
* ``options``, an `AiiDA`_ :py:class:`aiida.orm.nodes.data.Dict` instance, containing at least the keys ``resources``. More information about the options is available in the `AiiDA documentation`_.

Restarting calculations
^^^^^^^^^^^^^^^^^^^^^^^

The main difference between a ``VaspWorkChain`` and a  ``VaspCalculatio`` is that the former implements a basic logic of restarting failed or unfinished calculations.
The framework of :py:class:`~aiida.engine.processes.workchains.restart.BaseRestartWorkChain` is used with a set of predefined handlers to fix some (but not all) common pitfalls,
such as restarting an ionic relaxation that has run out of the wall time and electronic convergence issues. 

Once a calculation is finished, the ``CalculationNode`` is inspected by a series of :py:func:`aiida.engine.process.workchains.restart.process_handler`, 
which are executed in the order of descending priority.
Each handler may be tied to a specific list of ``exit_code`` that the calculation may have.
If any problems are found, and the restart can be performed, a ``ProcessHanlderReport`` would be returned and added to a list. 
If the ``break`` attribute of the report is set to ``True`` the handling process would be terminated.
Afterwards, the last report is inspected. If it has an none-zero ``exit_code`` the, then the workchain will be aborted with that ``exit_code`` returned, this corresponds to the case where the error cannot be corrected automatically.
Otherwise, it is assumed that calculation should be restarted with the revised inputs. 

The flow chart below illustrates how it works coupled with the emittion of the ``ProcessHanlderReport`` from the handlers:

.. image:: process-handler.png

For more information, please see the docstring of :py:class:`~aiida.engine.processes.workchains.restart.BaseRestartWorkChain`.

One should note that the hanlders included here are not intended to give a comprehensive coverage of all of possible errors from VASP, 
but instead we focus on improving the robustness by performing simple corrections that would be the right things to do in most times.

New handlers may be registered by adding the method to :py:class`~aiida_vasp.workchains.vasp.VaspWorkChain` with the ``process_handler`` decorator.
Alternatively, one can also extended the :py:class:`~aiida_vasp.workchains.vasp.VaspWorkChain` by sub-classing and add more handlers there.


.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/latest/
.. _Workchain: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/workflows.html#work-chains