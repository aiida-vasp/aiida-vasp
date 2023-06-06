.. _parameters:

Parameters
==========
Before describing how parameter passing works in this plugin it is worthwhile to restate that the design principle is that all higher lying workchains ultimately call the ``VaspWorkChain`` (:ref:`vasp_workchain`) which should handle `VASP`_ specific translations and setups in order to execute your problem with `VASP`_. At that point what we in general call parameters are fully converted to INCAR tags or flags in POSCAR, for instance in the case of selective dynamics.

.. note::
   In this documentation, there is the parameters, which is the general description of something you can adjust to get some specific behavior, or ``parameters`` which is
   a dedicated input parameter. It should be obvious from the context what is meant, and if in doubt, interpret the parameter as stated as the general meaning.

.. note::
   The higher lying workchains should in principle also be executable with e.g. Quantum Espresso or Abinit as long as one replaces the ``VaspWorkChain`` with a dedicated code specific translation and setup workchain for those codes. In order to facilitate this we designed what we call a ``ParameterMassager`` which translates certain parameters from AiiDA, or problem specific to `VASP`_ specific ones. This would also have to be written for the replacement code in case `VASP`_ is not to be used.

We now describe how parameters can be passed in the plugin. We separate between passing parameters directly to the ``VaspCalculation`` (:ref:`vasp_calculation`), the ``VaspWorkChain`` (or any workchain ultimately calling ``VaspWorkChain``). The latter being the recommended approach, unless you have very specific use-cases that warrants interacting with the ``VaspCalculation``.

When supplying parameters to ``VaspWorkChain``
----------------------------------------------
Any bundled workchain, :ref:`workchains` ultimately ends up calling ``VaspWorkChain`` and forwards parameters to it, where parameter composition is performed in terms dictated by `VASP`_ through the ``ParameterMassager``.
This is the recommended way of interacting and performing `VASP`_ calculations. For all parameters that go through the ``VaspWorkChain``, there are basically three approaches to supply and adjust parameters in the plugin:

Directly to VASP
^^^^^^^^^^^^^^^^
In case one needs to set certain parameters in ``INCAR`` that is unknown to the plugin, or one want to ensure VASP runs with some specific
tags in ``INCAR``, one can supply what we call *INCAR passthrough parameters* in ``parameters.incar``. The ``incar`` here is often referred to as a
namespace in ``parameters``. To assemble the ``parameters.incar`` using the ``ProcessBuilder`` do::

  builder = VaspWorkChain.get_builder()
  builder.parameters = Dict(dict={"incar": {"encut": 500, "prec": "accurate"}})
  ....
  submit(builder)

or through a plain ``dict`` definition converted to an AiiDA ``Dict``::

  inputs_dict = {"parameters":  Dict(dict={"incar": {"encut": 500, "prec": "accurate"}})}
  inputs_dict["structure"] = ...
  submit(VaspWorkChain, inputs_dict )

The relevant ``INCAR`` tags and values should be supplied as keys and values, respectively to the ``parameters.incar`` dictionary.
Remember that this is the only way to supply `VASP`_ parameters, or ``INCAR`` tags *directly* to `VASP`_ through the ``VaspWorkChain``, or
any more complex assembled workchain stack, eventually calling the ``VaspWorkChain``.

.. note::
   Notice that if in one of the workchains that composes the workchain stack, a parameter is detected that is incompatible with a successful execution of said workchain it is overridden,
   including any additional necessary parameters and a message to the user is provided in the logs giving details of the override.

.. note::
   Also note that for instance the degrees-of-freedom used for selective dynamics runs is to be considered as an input to `VASP`_ and this does for obvious reasons not fit into ``parameters.incar``.
   For these cases, a dedicated input node ``dynamics`` can be supplied to ``VaspWorkChain``.
   For this particular example, the key ``positions_dof`` contains the selective dynamics flags to be set in ``POSCAR``.
   Similarly one should utilize and define these namespaces accordingly when new functionality is introduced to `VASP`_ that for instance does not fit into ``INCAR``. Also, note that the parameter wording is not used for these inputs.

Workchain
^^^^^^^^^
When it comes to the workchains, one can supply parameters as direct inputs or contained in the ``parameters``. Take for instance the ``RelaxWorkChain``,
as an example. We can for the ``RelaxWorkChain`` either supply a direct input under the namespace ``relax``, i.e. ``relax.force_cutoff`` or supply it
in the ``parameters`` as ``parameters.relax.force_cutoff``.


The input it is also possible to supply *Workchain input parameter* in namespaces other than ``'incar'``, e.g. ``parameters.relax.force_cutoff``.
This would override any parameter that has been supplied (or set as default ) in e.g. ``parameters.incar.ediffg``.
This *Workchain input parameter* is translated to ``EDIFFG`` by the ``ParameterMassager`` which is called in the ``VaspWorkChain``.

Workflow
^^^^^^^^
From higher lying workchains like e.g. ``RelaxWorkChain``, the parameters pertaining to this workchain that would be relevant for `VASP`_ are passed to this workchain in a namespace called ``relax``, which is relevant if say the ``RelaxWorkChain`` is not the next callable workchain, but a few steps into the workflow. This generalizes to any workchain. This is again inspected in the relevant workchain and eventually forwarded to the next callable workchain in the workflow, eventually reaching ``VaspWorkChain``. As an example, such *Workflow input parameter* could be ``relax.force_cutoff`` is supplied, which controls the size of the force cutoff when performing relaxations and is used by the ``RelaxWorkChain`` or any other workchain along the way, which potentially need to modify this parameter. When leaving the ``RelaxWorkChain`` the relevant parameters with respect to relaxation should sit in ``parameters.relax``. The user can chose if they want to supply parameters in, for example ``relax`` or ``parameters.relax``.
Which strategy to follow depends on how explicit you want to be. For e.g. ``relax`` one typically define in the workchain ``spec.inputs`` section, ``relax.someparameter``.

A few general comments
^^^^^^^^^^^^^^^^^^^^^^
`VASP`_ specific parameters are only relevant for the first approach, while the two next approaches handle plugin specific parameters. For certain cases, like e.g. ``force_cutoff`` and negative ``EDIFFG`` there is correspondence, but this is not always the case.

The ``inherit_and_merge_parameters`` is called initially (before performing any workchain related work) in every workchain, which merges the workchain input nodes (``somenamespace.something``) with the content of ``parameters.somenamespace.something`` and performs the correct prioritization as described below. Notice that if an overlap is detected, the parameter in ``parameters.somenamespace.something`` is prioritized.

When supplying parameters to ``VaspCalculation``
-------------------------------------------------------------
The ``VaspCalculation`` expects the input of the `VASP`_ ``INCAR`` tags to be supplied using ``parameters``. E.g. the tags should be places accordingly, like ``parameters.icharg``, ``parameters.sigma`` etc. This is also what is provided from the ``parameters.incar`` output of the ``ParameterMassager``. In addition, the ``VaspCalculation`` accepts as a separate input node, the ``dynamics`` which can be supplied and contain e.g. the ``positions_dof`` flag which contains the selective dynamics flags used in ``POSCAR``.

Supported namespaces
--------------------
The supported namespaces are set using concatenation of the content of ``_BASE_NAMESPACES`` variable (currently containing ``['electronic', 'smearing', 'charge', 'dynamics', 'bands', 'relax', 'converge']``), any additional override namespace added by supplying the ``settings.additional_override_namespaces`` variable, which should be a list of strings to the ``VaspWorkChain`` and finally the override namespace ``incar`` directly related to `VASP`_ ``INCAR``.

How parameters are prioritized and set
--------------------------------------
There are basically three ways to supply parameters, some of them ultimately ending up as tags and values in the ``INCAR``, others are used to control the behavior of workchains or workflows composed of wotkchains. The prioritization of the different, but related ways are as follows:

1. *INCAR passthrough parameter*
2. *Workchain input parameter*
3. *Workflow input parameter*

When a conflict between these three ways to supply parameters occurs, the latter always overrides the formers.
For example, if you set (1) ``parameters.incar.ediffg = -1e-1``, (2) ``parameters.relax.force_cutoff = 1e-2``, and (3) ``relax.force_cutoff = Float(1e-3)`` for the ``RelaxWorkChain``, the parameter in (3) is finally chosen as `EDIFFG=-1e-3`.

The parameter massager
----------------------
The ``ParameterMassager`` translates parameters in the plugin to `VASP`_ specific ones, ensures that the prioritization is respected.
The input would be composed ``parameters`` containing elements from the workchain input nodes and the previously set ``parameters``, depending on the prioritization described above. The output would be a new ``parameters.incar`` which should only contain valid `VASP`_ flags, in addition to ``parameters.dynamics`` which should contain parameters to control the dynamics of the system (currently this only houses ``positions_dof`` to set the selective dynamics flags in ``POSCAR``. This is typically what is supplied to the ``VaspCalculation``.

Allowing custom `VASP`_ tags
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In case you for instance perform developments in the `VASP`_ code, sometimes it makes sense to add a new `VASP`_ tag. This can be supplied in ``settings.unsupported_parameters`` as dict with the following specifications::

  unsupported_parameters = {'my_unsupported_parameters': {
  'default': 1.0,
  'description': 'Some description',
  'type': float,
  'values': [1.0, 2.0]
  }
  builder.settings = Dict(dict={'unsupported_parameters': unsupported_parameters)

Alternatively, the validation can be turned off entirely by setting ``skip_parameters_validation`` to ``True`` under ``settings``, for example::

  builder.settings = Dict(dict={'skip_parameters_validation': True})

The above works for both ``VaspWorkChain`` and ``VaspCalculation``.
In the latter case, if any of ``skip_parameters_validation`` or ``unsupported_parameters`` are present in the ``settings`` input node, the validation is turned off completely.

.. _VASP: https://www.vasp.at
