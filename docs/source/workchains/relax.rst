.. _relax_workchain:

===============
Relax workchain
===============

Structure relaxation workchain, ``RelaxWorkChain``, which performs the regular duties of relaxing a structure. It has been designed such that calling workchains should try to use human readable parameters instead of the code dependent variables.

The workchain will try, if desired, to check the convergence of the structural relaxation, such as to try to ensure that the obtained result is not the result of the system being trapped in a local minima.

Running a single relaxation
----------------------------

An example containing necessary steps. Can be found in `run_relax`_ file in the examples folder. This includes creating sample input data and connecting it to the inputs of the workchain and finally how to submit the workchain The only steps not covered are the setup of AiiDA, which is documented in the `AiiDA docs`_

The example relaxes the ion positions only. Other degrees of freedom are available, see the detailed `input reference <relax_workchain/inputs>`

.. _run_relax: https://github.com/aiidateam/aiida-vasp/blob/develop/examples/run_relax.py
.. _AiiDA docs: https://aiida-core.readthedocs.io/en/stable/work/index.html

Writing workflows that include a relaxation step
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unless the workflow is solely concerned with providing a more robust or more specialized relaxation step, it is not recommended to derive from `vasp.relax`. Instead composing should be used; your workchain should set up the inputs and submit the relaxation workchain in a step function which will commit the results of the relaxation workchain to it's context. Further steps can then use the `relax.structure` result.

.. _relax_workchain/inputs:

Reference: `vasp.relax` inputs
------------------------------

The :py:class:`RelaxWorkChain<aiida_vasp.workchains.relax.RelaxWorkChain>` has several parameters that can be passed to it to control its behavior. Some are uniquely defined for the relaxation and others are exposed from the base `VaspWorkChain`.

.. note::
  Any parameter defined that would affect the relaxation, e.g. ``relax.positions`` will override any input passed via any other input, e.g passing ``ISIF`` manually in the ``parameters`` dictionary.


Basic inputs
^^^^^^^^^^^^

These are the set of basic parameters that need to be given to be able to perform the calculation. These inputs are required to perform the relaxation.

* ``structure``, type: :py:class:`aiida.orm.nodes.data.structure.StructureData` or :py:class:`aiida.orm.nodes.data.cif.CifData`. Describes the structure on which `VASP`_ is to be run.
* ``parameters``, type: :py:class:`aiida.orm.nodes.data.dict.Dict`. Dictionary with the parameters for the calculation. Please consult the documentation on how parameters are handled (:ref `parameters`) for details, particularly the section pertaining to the ``VaspWorkChain``.
* ``settings``, type: :py:class:`aiida.orm.nodes.data.dict.Dict`. Dictionary containing parameters not related to `VASP`_ itself, e.g. parser settings, selective dynamics, etc. **Optional**

Exposed from `vasp.vasp`
^^^^^^^^^^^^^^^^^^^^^^^^

The following inputs will be passed straight through to the :ref:`vasp_workchain`.

++++++++
Required
++++++++

These inputs are exposed from the base `VaspWorkChain`, and must be provided to perform the calculation.

* ``code``, type: :py:class:`aiida.orm.nodes.data.code.installed.InstalledCode`. Describes the VASP executable and holds a reference to the :py:class:`Computer<aiida.orm.computers.Computer>` instance on which it lives.
* ``kpoints``, type: :py:class:`aiida.orm.nodes.data.array.kpoints.KpointsData`. The kpoints mesh or path for the calculation.
* ``potential_family``, type: :py:class:`aiida.orm.nodes.data.str.Str`. The name given to a set of uploaded POTCAR files.
* ``potential_mapping``, type: :py:class:`aiida.orm.nodes.data.dict.Dict`. Dictionary containing an entry for at least every kind name in the ``structure`` input with the full name of the POTCAR from the ``potential_family``. Example: ``{'In1': 'In_d', 'In2': 'In_h'}``.
* ``options``, type: :py:class:`aiida.orm.nodes.data.dict.Dict`. Dictionary containing at least the keys ``resources``. More information about the options is available in the `AiiDA documentation`_.

++++++
Extras
++++++

These inputs do not need to be provided and have a set of defaults.

* ``max_iterations``, type: :py:class:`aiida.orm.nodes.data.int.Int`, default: 5. How many iterations the restart will be attempted before resulting in failure. -> `max_iterations` on `vasp.vasp`
* ``clean_workdir``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: True. Whether or not the remote folder of the calculation will be deleted after the end of the calculation. -> `clean_workdir` on `vasp.vasp`

Relaxation control
^^^^^^^^^^^^^^^^^^

.. _EDIFFG: https://www.vasp.at/wiki/index.php/EDIFFG
.. _EDIFF: https://www.vasp.at/wiki/index.php/EDIFF

These inputs control global parameters about the relaxation.

All of these inputs are optional

* ``relax.perform``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: False. Whether or not to perform relaxations
* ``relax.steps``, type: :py:class:`aiida.orm.nodes.data.int.Int`, default: 60. The number of ionic positions updates to perform.
* ``relax.keep_magnetization``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: True. Whether or not to keep the magnetization from the previous relaxation run.
* ``relax.algo``, type: :py:class:`aiida.orm.nodes.data.str.Str`, default: cg. The type of algorithm that will be used for the ionic relaxation.
* ``relax.energy_cutoff``, type: :py:class:`aiida.orm.nodes.data.float.Float`. The cutoff that determines when the relaxation procedure is stopped. In this case it stops when the total energy between two ionic steps is less than the supplied value. If not provided whatever default value `VASP`_ has for `EDIFF`_.
* ``relax.force_cutoff``, type: :py:class:`aiida.orm.nodes.data.float.Float`. The cutoff that determines when the relaxation procedure is stopped. In this case it stops when all forces are smaller than than the supplied value. If not provided whatever default value `VASP`_ has for `EDIFFG`_.
* ``relax.perform_static``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: True. Whether or not to perform a static calculation after the relaxation.

Degrees of freedom
^^^^^^^^^^^^^^^^^^

Through its input parameters, `VASP`_ can be configured to utilize three degrees of freedom for relaxations: ion positions, cell volume and cell shape. Some, but not all combinations are allowed, read more about that in the `official VASP wiki - ISIF tag page`_. Other possibilities are also doable, but not covered here and typically demands a dedicated VASP version with hard coded changes to the source code.

`vasp.relax` allows to switch each degree of freedom on / off independently, setting the ``ISIF`` and ``IBRION`` tags accordingly. Each of these inputs is optional and by default only the ion positions are relaxed.

All of these inputs are optional

* ``relax.positions``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: True. If is True, perform relaxations of the atomic positions.
* ``relax.shape``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: False. If is True, perform relaxation of the cell shape.
* ``relax.volume``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: False. If is True, perform relaxation of the cell volume.

.. _official VASP wiki - ISIF tag page: https://cms.mpi.univie.ac.at/wiki/index.php/ISIF

Convergence
^^^^^^^^^^^

The `vasp.relax` workchain is able to check for convergence on any of the available degrees of freedom by running a fresh relaxation from the output structure of the previous run. This can sometimes lead to further relaxation, if the previous run got stuck in a local charge density minimum. The new calculation starts from scratch with a randomized charge density, but with the last obtained positions. This is done iteratively until the target property does not change more than a given tolerance. Currently, external check on force, stress and energy is not implemented, but will be available as an option in the future.

Keep in mind there is no guarantee that the new run will overcome the barriers of a local minimum. More in-depth workchains could be developed to do that, by deriving from this workchain or using it as a building block. This feature is switched off by default.

All of these inputs are optional

* ``relax.convergence_on``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: False. Whether or not to check or run additional relaxations.
* ``relax.convergence_absolute``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: False. Whether or not to converge the relaxation with respect to the previous run
  - False: relative tolerances are used (relative convergence)
  - True: absolute tolerances are used (native VASP units)
* ``relax.convergence_max_iterations``, type: :py:class:`aiida.orm.nodes.data.int.Int`, default: 5. Maximum number of relaxation runs.
* ``relax.convergence_shape_lengths``, type: :py:class:`aiida.orm.nodes.data.float.Float`, default: 0.1. Maximum percentage change of the L2 norm for the unitcell vectors from the previous relaxation.
* ``relax.convergence_shape_angles``, type: :py:class:`aiida.orm.nodes.data.float.Float`, default: 0.1. Maximum percentage change of the unitcell angles from the previous relaxation.
* ``relax.convergence_volume``, type: :py:class:`aiida.orm.nodes.data.float.Float`, default: 0.01. Maximum percentage change of the unitcell volume from the previous relaxation.
* ``relax.convergence_positions``, type: :py:class:`aiida.orm.nodes.data.float.Float`, default: 0.01. Maximum percentage change of the displacement (L2 norm) of the positions from the previous relaxation.

Reference: `vasp.relax` outputs
-------------------------------

The following output nodes are created upon successful completion:

* ``misc``, type: :py:class:`aiida.orm.nodes.data.dict.Dict`. Dictionary containing the output parameters containing smaller quantities that do not depend on system size.
* ``relax.structure``, type: :py:class:`aiida.orm.nodes.data.structure.StructureData`. The output structure after relaxation (if it was performed).

Depending on the passed inputs to the workchain several outputs might be exposed according to what was defined in the :ref:`vasp_workchain_outputs`.

.. _VASP: https://www.vasp.at
.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/latest/
