.. _parameters:

Parameters
==========
Before describing how parameter passing works in this plugin it is worthwhile to restate that the design principle is that all higher lying workchains ultimately call the ``VaspWorkChain`` (:ref:`vasp_workchain`) which should handle `VASP`_ specific translations and setups in order to execute your problem with `VASP`_. The higher lying workchains should in principle also be executable with e.g. Quantum Espresso or Abinit as long as one replaces the ``VaspWorkChain`` with a dedicated code specific translation and setup workchain for those codes. In order to facilitate this we designed what we call a ``ParameterMassager`` which translates certain parameters from AiiDA to `VASP`_ specific ones. This would also have to be written for the replacement code.

We now describe how parameters can be passed in the plugin. We separate between passing parameters directly to the ``VaspCalculation`` (:ref:`vasp_calculation`), the ``VaspWorkChain`` (or any workchain ultimately calling ``VaspWorkChain``).

When supplying parameters to ``VaspWorkChain``
----------------------------------------------
Also valid for any workchain (:ref:`workchains`) ultimately ending up calling ``VaspWorkChain`` and forwarding parameters to it. This is the recommended way of interacting and performing `VASP`_ calculations. Here, one can supply `VASP`_ input parameters (``INCAR`` tags) in ``inputs.parameters.incar``. Any key value combination put into this ``incar`` namespace bypasses the ``ParameterMassager`` (except the check if it is valid `VASP`_ tags). Remember that this is the only way to supply `VASP`_ parameters, or ``INCAR`` tags directly to `VASP`_ through the workchain stack. In the ``ParameterMassager`` other namespaces can also be supplied. For instance from higher lying workchains like the ``RelaxWorkChain``, the parameters pertaining to this workchain that would be relevant for `VASP`_ are passed to the ``VaspWorkChain`` in a namespace called ``relax`` and in there parameters like ``force_cutoff`` is supplied which controls the size of the force cutoff when performing relaxations. This parameter is translated to a suitable ``EDIFFG`` by the ``ParameterMassager``. To the workchains, e.g. the ``RelaxWorkChain`` it is also possible to supply these namespaces in the ``inputs.parameters``, e.g. ``inputs.parameters.relax.force_cutoff`` and this would override any parameter that has been supplied (or set as default ) in ``input.relax.force_cutoff`` (an input node of the ``RelaxWorkChain``). In a workchain, the ``inherit_and_merge_parameters`` is called, which before the next workchain is called merges the workchain input nodes (``inputs.somenamespace.something``) with the content of ``inputs.parameters.somenamespace.something`` and performs the correct prioritization as described below.

When supplying parameters directly to the ``VaspCalculation``
-------------------------------------------------------------
The ``VaspCalculation`` expects the input of the `VASP`_ ``INCAR`` tags to be supplied using ``inputs.parameters``. E.g. the tags should be places accordingly, like ``inputs.parameters.icharg``, ``inputs.parameters.sigma`` etc. This is also what is provided from the ``parameters.incar`` output of the ``ParameterMassager``. In addition, the ``VaspCalculation`` accepts as a separate input node, the ``inputs.dynamics`` which can be supplied and contain e.g. the ``positions_dof`` flag which contains the selective dynamics flags used in ``POSCAR``.

Supported namespaces
--------------------
The supported namespaces are set using concatenation of the content of ``_BASE_NAMESPACES`` variable (currently containing ``['electronic', 'smearing', 'charge', 'dynamics', 'bands', 'relax', 'converge']``), any additional override namespace added by supplying the ``inputs.settings.additional_override_namespaces`` variable, which should be a list of strings to the ``VaspWorkChain`` and finally the override namespace ``incar``. 

How parameters are prioritized and set
--------------------------------------
Since there are basically three ways to supply parameters (for the most general ones ultimately ending up as tags in the ``INCAR``),
the following prioritization is performed:

1. Parameters supplied in the ovveride namespace ``inputs.parameters.incar`` are always respected as long as it is a valid `VASP`_ tag.
2. Parameters supplied in ``inputs.parameters``, e.g. ``inputs.parameters.relax.force_cutoff`` (for the ``RelaxWorkChain``).
3. Parameters supplied as workchain input nodes, e.g. ``inputs.relax.force_cutoff`` (for the ``RelaxWorkChain``).

The parameter massager
----------------------
The ``ParameterMassager`` translates parameters in the plugin to `VASP`_ specific ones, ensures that the prioritization is respected.
The input would be composed ``inputs.parameters`` containing elements from the workchain input nodes and the previously set ``inputs.parameters``, depending on the prioritization described above. The output would be a new ``parameters.incar`` which should only contain valid `VASP`_ flags, in addition to ``parameters.dynamics`` which should contain parameters to control the dynamics of the system (currently this only houses ``positions_dof`` to set the selective dynamics flags in ``POSCAR``. This is typically what is supplied to the ``VaspCalculation``.

Allowing custom `VASP`_ tags
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In case you for instance perform developments in the `VASP`_ code, sometimes it makes sense to add a new `VASP`_ tag. This can be supplied in ``settings.inputs.unsupported_parameters`` as dict with the following specifications::

  unsupported_parameters = {'my_unsupported_parameters': {
  'default': 1.0,
  'description': 'Some description',
  'type': float,
  'values': [1.0, 2.0]
  }

.. _VASP: https://www.vasp.at
