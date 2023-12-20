.. _converge_workchain:

==================
Converge workchain
==================

The :py:class:`ConvergeWorkChain<aiida_vasp.workchains.converge.ConvergeWorkChain>` is intended to be used to determine the appropriate energy cutoff and reciprocal space mesh as to converge a certain convergence quantity (``cutoff_type``), e.g. the ground state total energy. One can control whether only the energy, reciprocal mesh or both parameters are used.

When performing the energy cutoff convergence if no k-point mesh is given the energy cutoff convergence tests will be done with an auto-generated k-point mesh. The distance between the points is either user provided or default values are used. Similarly, for the reciprocal mesh, if a provided energy cutoff has been given that will be used, otherwise if a cutoff value has been found from the energy convergence that would be used instead.

At the end of each convergence test the value for the energy cutoff and/or k-points grid will be taken and stored.

For example, lets say that the user provides the k-points but no energy cutoff. The convergence test will run a series of calculations (serially) where the value of the energy cutoff is varied. After each calculation the workflow will try to see if the desired convergence criterion was reached. If it was the workflow will perform a last calculation with the recommended value for the energy cutoff and the provided k-points, and will store these values for further usage in other calculations.

The workflow works in the same way for the k-points, if the cutoff energy is given the same procedure will be done, but varying the number of k-points.

The situation is somewhat more complex when both variables are allowed to vary, since first a convergence in energy will be performed with a k-point mesh given by `converge.k_spacing`, which is deemed to be "good enough". After the cutoff energy recommended value has been found, that value will be set and the k-points are varied instead. Once convergence has been reached both values will be set and the final calculation will be performed and the recommended values stored for further use.

.. note::
  The details of the :py:class:`ConvergeWorkChain<aiida_vasp.workchains.converge.ConvergeWorkChain>` can be also seen in the ``run_converge.py`` example.


Reference: `vasp.converge` inputs
---------------------------------

These are the set of basic parameters for the base level convergence of the energy cutoff and reciprocal space mesh.

* ``parameters``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary with the parameters for the calculation. Please consult the documentation on how parameters are handled (:ref `parameters`) for details, particularly the section pertaining to the ``VaspWorkChain``.
* ``structure``, type: :py:class:`StructureData<aiida.orm.nodes.data.structure.StructureData>` or :py:class:`CifData<aiida.orm.nodes.data.cif.CifData>`. Describes the structure on which `VASP`_ is to be run.
* ``settings``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary containing parameters not related to `VASP`_ itself, e.g. parser settings, selective dynamics, etc. **Optional**
* ``kpoints``, type :py:class:`KpointsData<aiida.orm.nodes.data.array.kpoints.KpointsData>`. If given no convergence will be performed in the reciprocal mesh. Instead this mesh will be used while the energy cutoff is varied. **Optional**

Convergence parameters
^^^^^^^^^^^^^^^^^^^^^^

These parameters dictate the behavior of the workchain during the convergence tests. These can be used to control which kind of convergence test will be done, how strict, etc.

All these parameters are optional.

* ``converge.pwcutoff``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`. The plane-wave cutoff to be used during convergence tests in electron volts.
* ``converge.kgrid``, type: :py:class:`ArrayData<aiida.orm.nodes.data.array.array.ArrayData>`. The k-point grid to be used during convergence tests.
* ``converge.pwcutoff_start``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`, default: 200. The initial value of the plane-wave cutoff in electron volts for the convergence tests.
* ``converge.pwcutoff_step``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`, default: 50. The plane-wave cutoff step (increment) in electron volts.
* ``converge.pwcutoff_samples``, type: :py:class:`Int<aiida.orm.nodes.data.int.Int>`, default: 10. The number of plane-wave cutoff samples.
* ``converge.k_dense``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`, default: 0.07. The target k-point stepping at the densest grid in inverse Ångströms.
* ``converge.k_coarse``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`, default: 0.35. The target k-point stepping at the coarsest grid in inverse Ångströms.
* ``converge.k_spacing``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`, default: 0.1. The default k-point spacing in inverse Ångströms. This value will be used to perform the energy cutoff convergence test if no k-points are explicitly given.
* ``converge.k_samples``, type: :py:class:`Int<aiida.orm.nodes.data.int.Int>`, default: 10. The number of k-point samples.
* ``converge.cutoff_type``, type: :py:class:`Str<aiida.orm.nodes.data.str.Str>`, default: energy. The cutoff_type to check convergence against. Currently the following options are accepted: energy, forces, gap and vbm (not yet currently supported).
* ``converge.cutoff_value``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`, default: 0.1. The cutoff value to be used when considering absolute differences. When the difference between two convergence calculations are within this value for ``cutoff_type``, then it is considered converged.
* ``converge.cutoff_value_r``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`, default: 0.01. The relative cutoff value to be used. When the difference between two convergence calculations are within this value for ``cutoff_type``, then it is considered converged. However, in this case the cutoff value is the difference between `cutoff_type` for the input structure and an atomic displacement or a compression of the unitcell.
* ``converge.compress``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: False. If True, a convergence test of the compressed structure is also performed. The difference of the ``cutoff_type`` values for each calculations are evaluated and when the difference between these are less than ``cutoff_value_r``, the calculation is considered converged. The largest plane-wave cutoff and densest k-point grid are used.
* ``converge.displace``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: False. If True, a convergence test of the displaced structure is also performed. The difference of the ``cutoff_type`` values for each calculations are evaluated and when the difference between these are less than ``cutoff_value_r``, the calculation is considered converged. The largest plane-wave cutoff and densest k-point grid are used.
* ``converge.displacement_vector``, type: :py:class:`ArrayData<aiida.orm.nodes.data.array.array.ArrayData>`, default: [1.0, 1.0, 1.0]. The displacement unit vector for the displacement test. Sets the direction of displacement.
* ``converge.displacement_atom``, type: :py:class:`Int<aiida.orm.nodes.data.int.Int>`, default: 1. Which atom to displace? Index starts from 1 and follows the sequence for the sites in the Aiida ``structure`` object.
* ``converge.volume_change``, type: :py:class:`ArrayData<aiida.orm.nodes.data.array.array.ArrayData>`, default: [1.05, 1.05, 1.05]. The volume change in direct coordinates for each lattice vector.
* ``converge.relax``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: False. Whether or not to perform a relaxation for each convergence test.
* ``converge.total_energy_type``, type: :py:class:`Str<aiida.orm.nodes.data.str.Str>`, default: energy_extrapolated. The energy type that is used when ``cutoff_type`` is set to `energy`. Consult the options available in the parser for the current version.
* ``converge.testing``,type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: False. Whether or not the calculation is a test. Mostly used for debugging and CI/CD.

Exposed from `vasp.vasp`
^^^^^^^^^^^^^^^^^^^^^^^^

The following inputs will be passed straight through to the :ref:`vasp_workchain`.

Required
""""""""

These inputs are exposed from the base `VaspWorkChain`, and must be provided to perform the calculation.

* ``code``, type: :py:class:`InstalledCode<aiida.orm.nodes.data.code.installed.InstalledCode>`. Describes the VASP executable and holds a reference to the :py:class:`Computer<aiida.orm.computers.Computer>` instance on which it lives.
* ``potential_family``, type: :py:class:`Str<aiida.orm.nodes.data.str.Str>`. The name given to a set of uploaded POTCAR files.
* ``potential_mapping``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary containing an entry for at least every kind name in the ``structure`` input with the full name of the POTCAR from the ``potential_family``. Example: ``{'In1': 'In_d', 'In2': 'In_h'}``.
* ``options``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary containing at least the keys ``resources``. More information about the options is available in the `AiiDA documentation`_.

Extras
""""""

These inputs do not need to be provided and have a set of defaults.

* ``max_iterations``, type: :py:class:`Int<aiida.orm.nodes.data.int.Int>`, default: 5. How many iterations the restart will be attempted before resulting in failure. -> `max_iterations` on `vasp.vasp`
* ``clean_workdir``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: True. Whether or not the remote folder of the calculation will be deleted after the end of the calculation. -> `clean_workdir` on `vasp.vasp`
* ``settings``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary containing parameters not related to `VASP`_ itself, e.g. parser settings, selective dynamics, etc.

Exposed from `vasp.relax`
^^^^^^^^^^^^^^^^^^^^^^^^^

.. _EDIFFG: https://www.vasp.at/wiki/index.php/EDIFFG
.. _EDIFF: https://www.vasp.at/wiki/index.php/EDIFF
.. _official VASP wiki - ISIF tag page: https://cms.mpi.univie.ac.at/wiki/index.php/ISIF

These inputs control global parameters about the relaxation. These are passed to the underlying `RelaxWorkChain` which is called during each step of the `ConvergeWorkChain`. Whether or not an actual relaxation if performed depends on the value of ``converge.relax`` .

All of these inputs are optional

* ``relax.perform``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: False. Whether or not to perform relaxations
* ``relax.positions``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: True. If is True, perform relaxations of the atomic positions.
* ``relax.shape``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: False. If is True, perform relaxation of the cell shape.
* ``relax.volume``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: False. If is True, perform relaxation of the cell volume.
* ``relax.steps``, type: :py:class:`Int<aiida.orm.nodes.data.int.Int>`, default: 60. The number of ionic positions updates to perform.
* ``relax.keep_magnetization``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: True. Whether or not to keep the magnetization from the previous relaxation run.
* ``relax.algo``, type: :py:class:`Str<aiida.orm.nodes.data.str.Str>`, default: cg. The type of algorithm that will be used for the ionic relaxation.
* ``relax.energy_cutoff``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`. The cutoff that determines when the relaxation procedure is stopped. In this case it stops when the total energy between two ionic steps is less than the supplied value. If not provided whatever default value `VASP`_ has for `EDIFF`_.
* ``relax.force_cutoff``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`. The cutoff that determines when the relaxation procedure is stopped. In this case it stops when all forces are smaller than than the supplied value. If not provided whatever default value `VASP`_ has for `EDIFFG`_.
* ``relax.convergence_on``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: False. Whether or not to check or run additional relaxations.
* ``relax.convergence_absolute``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: False. Whether or not to converge the relaxation with respect to the previous run
  - False: relative tolerances are used (relative convergence)
  - True: absolute tolerances are used (native VASP units)
* ``relax.convergence_max_iterations``, type: :py:class:`Int<aiida.orm.nodes.data.int.Int>`, default: 5. Maximum number of relaxation runs.
* ``relax.convergence_shape_lengths``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`, default: 0.1. Maximum percentage change of the L2 norm for the unitcell vectors from the previous relaxation.
* ``relax.convergence_shape_angles``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`, default: 0.1. Maximum percentage change of the unitcell angles from the previous relaxation.
* ``relax.convergence_volume``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`, default: 0.01. Maximum percentage change of the unitcell volume from the previous relaxation.
* ``relax.convergence_positions``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`, default: 0.01. Maximum percentage change of the displacement (L2 norm) of the positions from the previous relaxation.
* ``perform_static``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: True. Whether or not to perform a static calculation after the relaxation.

Reference: `vasp.converge` outputs
----------------------------------

The following output nodes are created upon successful completion:

* ``misc``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary containing the output parameters containing smaller quantities that do not depend on system size.

Extras
^^^^^^

These outputs might be present depending on the type of calculation performed, i.e. which kind of convergence, if any, was performed.

* ``converge.data``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary containing the value of the convergence criterion parameter for each variation of the convergence parameters (energy and/or k-points).
* ``converge.pwcutoff_recommended``, type: :py:class:`Float<aiida.orm.nodes.data.float.Float>`. Recommended value for the energy cutoff.
* ``converge.kpoints_recommended``, type: :py:class:`KpointsData<aiida.orm.nodes.data.array.kpoints.KpointsData>`. Recommended value for the k-points mesh.

.. _VASP: https://www.vasp.at
.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/latest/
