.. _converge_workchain:

==================
Converge workchain
==================

The ``ConvergeWorkChain`` is intended to be used to determine the appropriate energy cutoff and reciprocal space mesh as to converge the ground state total energy. One can control whether only the energy, reciprocal mesh or both parameters are used.

When performing the energy cutoff convergence if no k-point mesh is given the energy cutoff convergence tests will be done with an auto-generated k-point mesh. The distance between the points is either user provided or default values are used. Similarly, for the reciprocal mesh, if a provided energy cutoff has been given that will be used, otherwise if a cutoff value has been found from the energy convergence that would be used instead.

.. note::
   The details of the ``ConvergeWorkChain`` can be also seen in the ``run_converge.py`` example.


Reference: `vasp.converge` inputs
---------------------------------

These are the set of basic parameters for the base level convergence of the energy cutoff and reciprocal space mesh.

* ``parameters``, type: :py:class:`aiida.orm.nodes.data.Dict`. Dictionary with the parameters for the calculation. Please consult the documentation on how parameters are handled (:ref `parameters`) for details, particularly the section pertaining to the ``VaspWorkChain``.
* ``structure``, type: :py:class:`aiida.orm.nodes.data.StructureData` or :py:class:`aiida.orm.nodes.data.CifData`. Describes the structure on which `VASP`_ is to be run.


Convergence parameters
^^^^^^^^^^^^^^^^^^^^^^

These parameters dictate the behaviour of the workchain during the convergence tests. These can be used to control which kind of convergence test will be done, how strict, etc.

All these 

* ``converge.pwcutoff``, type: :py:class:`aiida.orm.nodes.data.Float`. The plane-wave cutoff to be used during convergence tests in electron volts.
* ``converge.kgrid``, type: :py:class:`aiida.orm.nodes.data.ArrayData`. The k-point grid to be used during convergence tests.
* ``converge.pwcutoff_start``, type: :py:class:`aiida.orm.nodes.data.Float`, default: 200. The initial value of the plane-wave cutoff in electron volts for the convergence tests.
* ``converge.pwcutoff_step``, type: :py:class:`aiida.orm.nodes.data.Float`, default: 50. The plane-wave cutoff step (increment) in electron volts.
* ``converge.pwcutoff_samples``, type: :py:class:`aiida.orm.nodes.data.Int`, default: 10. The number of plane-wave cutoff samples.
* ``converge.k_dense``, type: :py:class:`aiida.orm.nodes.data.Float`, default: 0.07. The target k-point stepping at the densest grid in inverse AA.
* ``converge.k_course``, type: :py:class:`aiida.orm.nodes.data.Float`, default: 0.35. The target k-point stepping at the coarsest grid in inverse AA.
* ``converge.k_spacing``, type: :py:class:`aiida.orm.nodes.data.Float`, default: 0.1. The default k-point spacing in inverse AA.
* ``converge.k_samples``, type: :py:class:`aiida.orm.nodes.data.Int`, default: 10. The number of k-point samples.

Exposed from `vasp.vasp`
^^^^^^^^^^^^^^^^^^^^^^^^

The following inputs will be passed straight through to the :ref:`vasp_workchain`.

++++++++
Required
++++++++

These inputs are exposed from the base `VaspWorkChain`, and must be provided to perform the calculation.

* ``code``, type: :py:class:`aiida.orm.nodes.data.Code`. Describes the VASP executable and holds a reference to the ``Computer`` instance on which it lives.
* ``potential_family``, type: :py:class:`aiida.orm.nodes.data.Str`. The name given to a set of uploaded POTCAR files.
* ``potential_mapping``, type: :py:class:`aiida.orm.nodes.data.Dict`. Dictionary containing an entry for at least every kind name in the ``structure`` input with the full name of the POTCAR from the ``potential_family``. Example: ``{'In1': 'In_d', 'In2': 'In_h'}``.
* ``options``, type: :py:class:`aiida.orm.nodes.data.Dict`. Dictionary containing at least the keys ``resources``. More information about the options is available in the `AiiDA documentation`_.

++++++
Extras
++++++

These inputs do not need to be provided and have a set of defaults.

* ``max_iterations``, type: :py:class:`aiida.orm.nodes.data.Int`, default: 5. How many iterations the restart will be attempted before resulting in failure. -> `max_iterations` on `vasp.vasp`
* ``clean_workdir``, type: :py:class:`aiida.orm.nodes.data.Bool`, default: True. Whether or not the remote folder of the calculation will be deleted after the end of the calculation. -> `clean_workdir` on `vasp.vasp`
* ``settings``, type: :py:class:`aiida.orm.nodes.data.Dict`. Dictionary containing parameters not related to `VASP`_ itself, e.g. parser settings, selective dynamics, etc.

Exposed from `vasp.relax`
^^^^^^^^^^^^^^^^^^^^^^^^^

.. _EDIFFG: https://www.vasp.at/wiki/index.php/EDIFFG
.. _EDIFF: https://www.vasp.at/wiki/index.php/EDIFF
.. _official VASP wiki - ISIF tag page: https://cms.mpi.univie.ac.at/wiki/index.php/ISIF

These inputs control global parameters about the relaxation. These are passed to the underlying `RelaxWorkChain` which is called during each step of the `ConvergeWorkChain`.

All of these inputs are optional

* ``relax.perform``, type: :py:class:`aiida.orm.nodes.data.Bool`, default: False. Whether or not to perform relaxations
* ``relax.positions``, type: :py:class:`aiida.orm.nodes.data.Bool`, default: True. If is True, perform relaxations of the atomic positions.
* ``relax.shape``, type: :py:class:`aiida.orm.nodes.data.Bool`, default: False. If is True, perform relaxation of the cell shape.
* ``relax.volume``, type: :py:class:`aiida.orm.nodes.data.Bool`, default: False. If is True, perform relaxation of the cell volume.
* ``relax.steps``, type: :py:class:`aiida.orm.nodes.data.Int`, default: 60. The number of ionic positions updates to perform.
* ``relax.keep_magnetization``, type: :py:class:`aiida.orm.nodes.data.Bool`, default: True. Whether or not to keep the magnetization from the previous relaxation run.
* ``relax.algo``, type: :py:class:`aiida.orm.nodes.data.Str`, default: cg. The type of algorithm that will be used for the ionic relaxation.
* ``relax.energy_cutoff``, type: :py:class:`aiida.orm.nodes.data.Float`. The cutoff that determines when the relaxation procedure is stopped. In this case it stops when the total energy between two ionic steps is less than the supplied value. If not provided whatever default value `VASP`_ has for `EDIFF`_.
* ``relax.force_cutoff``, type: :py:class:`aiida.orm.nodes.data.Float`. The cutoff that determines when the relaxation procedure is stopped. In this case it stops when all forces are smaller than than the supplied value. If not provided whatever default value `VASP`_ has for `EDIFFG`_.
* ``relax.convergence_on``, type: :py:class:`aiida.orm.nodes.data.Bool`, default: False. Whether or not to check or run additional relaxations.
* ``relax.convergence_absolute``, type: :py:class:`aiida.orm.nodes.data.Bool`, default: False. Whether or not to converge the relaxation with respect to the previous run
  - False: relative tolerances are used (relative convergence)
  - True: absolute tolerances are used (native VASP units)
* ``relax.convergence_max_iterations``, type: :py:class:`aiida.orm.nodes.data.Int`, default: 5. Maximum number of relaxation runs.
* ``relax.convergence_shape_lengths``, type: :py:class:`aiida.orm.nodes.data.Float`, default: 0.1. Maximum percentage change of the L2 norm for the unitcell vectors from the previous relaxation.
* ``relax.convergence_shape_angles``, type: :py:class:`aiida.orm.nodes.data.Float`, default: 0.1. Maximum percentage change of the unitcell angles from the previous relaxation.
* ``relax.convergence_volume``, type: :py:class:`aiida.orm.nodes.data.Float`, default: 0.01. Maximum percentage change of the unitcell volume from the previous relaxation.
* ``relax.convergence_positions``, type: :py:class:`aiida.orm.nodes.data.Float`, default: 0.01. Maximum percentage change of the displacement (L2 norm) of the positions from the previous relaxation.
* ``perform_static``, type: :py:class:`aiida.orm.nodes.data.Bool`, default: True. Whether or not to perform a static calculation after the relaxation.