.. _bands_workchain:

===============
Bands workchain
===============

The ``BandsWorkChain`` deals with the extraction of the band structure of a material using `seekpath`_ as a pre-processor to get the k-points path.

This workchain will take the results from a previous ground state calculation and use them to determine the band structure itself.

.. note::
   The details of the :py:class:`BandsWorkChain<aiida_vasp.workchains.bands.BandsWorkChain>` can be also seen in the ``run_bands.py`` example.


Reference: `vasp.bands` inputs
------------------------------

Basic inputs
^^^^^^^^^^^^

These are the set of basic parameters for the base level determination of the band structure.

++++++++
Required
++++++++

The calculation of the band structure assumes that a previous calculation where the ground state has already been performed. The results of that calculation must be stored in a :py:class:`aiida.orm.nodes.data.remote.base.RemoteData` folder resulting from a previous `aiida-vasp` run.

* ``restart_folder``, type: :py:class:`aiida.orm.nodes.data.remote.base.RemoteData`. The folder to restart in, which contains the outputs from the previously preformed calculation to extract the charge density.

++++++
Extras
++++++

These inputs do not need to be provided for the band structure calculation.

* ``parameters``, type: :py:class:`aiida.orm.nodes.data.dict.Dict`. Dictionary with the parameters for the calculation. Please consult the documentation on how parameters are handled (:ref `parameters`) for details, particularly the section pertaining to the ``VaspWorkChain``.
* ``settings``, type: :py:class:`aiida.orm.nodes.data.dict.Dict`. Dictionary containing parameters not related to `VASP`_ itself, e.g. parser settings, selective dynamics, etc.

Exposed from `vasp.vasp`
^^^^^^^^^^^^^^^^^^^^^^^^

The following inputs will be passed straight through to the :ref:`vasp_workchain`.

++++++++
Required
++++++++

These inputs are exposed from the base `VaspWorkChain`, and must be provided to perform the calculation.

* ``structure``, type: :py:class:`aiida.orm.nodes.data.structure.StructureData` or :py:class:`aiida.orm.nodes.data.CifData`. Describes the structure on which `VASP`_ is to be run.
* ``code``, type: :py:class:`aiida.orm.nodes.data.Code`. Describes the VASP executable and holds a reference to the ``Computer`` instance on which it lives.
* ``potential_family``, type: :py:class:`aiida.orm.nodes.data.str.Str`. The name given to a set of uploaded POTCAR files.
* ``potential_mapping``, type: :py:class:`aiida.orm.nodes.data.dict.Dict`. Dictionary containing an entry for at least every kind name in the ``structure`` input with the full name of the POTCAR from the ``potential_family``. Example: ``{'In1': 'In_d', 'In2': 'In_h'}``.
* ``options``, type: :py:class:`aiida.orm.nodes.data.dict.Dict`. Dictionary containing at least the keys ``resources``. More information about the options is available in the `AiiDA documentation`_.

++++++
Extras
++++++

These inputs do not need to be provided and have a set of defaults.

* ``max_iterations``, type: :py:class:`aiida.orm.nodes.data.int.Int`, default: 5. How many iterations the restart will be attempted before resulting in failure. -> `max_iterations` on `vasp.vasp`
* ``clean_workdir``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: True. Whether or not the remote folder of the calculation will be deleted after the end of the calculation. -> `clean_workdir` on `vasp.vasp`

Smearing
^^^^^^^^

These parameters control the smearing of the charge density when determining the band structure.

* ``smearing.gaussian``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: True. Whether or not gaussian smearing would be used in the band structure determination. If it is not set the smearing would be set to Fermi smearing.
* ``smearing.sigma``, type: :py:class:`aiida.orm.nodes.data.float.Float`, default: 0.05. Magnitude of the smearing applied to the band structure determination, in eV.

Bands specific information
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _LORBIT: https://www.vasp.at/wiki/index.php/LORBIT

The following inputs will allow the user to control how the band structure is determined in greater detail. None of these inputs are required and all have pre-defined default values.

* ``bands.kpoints_distance``, type: :py:class:`aiida.orm.nodes.data.float.Float`, default: 0.05. The distance between each k-point along each high-symmetry line.
* ``bands.decompose_bands``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: False. Whether or not the bands will be decomposed per atom.
* ``bands.decompose_wave``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: False. Whether to decompose the wave function when determining the band structure.
* ``bands.lm``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: False. Whether or not to decompose the wave function into l- and m- states.
* ``bands.phase``, type: :py:class:`aiida.orm.nodes.data.bool.Bool`, default: False. Whether or not to further decompose the l- and m- projections into phases.
* ``bands.wigner_seitz_radius``, type: :py:class:`aiida.orm.nodes.data.list.List`, default: `list[False]`. The Wigner-Seitz radius for each atom type in ångstroms as a list. If set, the internal projectors are not utilized.

.. note::
   The parameters dealing with the decomposition of the wave function, ``bands.decompose_bands``, ``bands.decompose_wave``, etc. will be used to determine the value for `LORBIT`_ needed to fulfill the desired decompositions.

   These bands specific values will override any value passed via the ``parameters``, e.g. `LORBIT`_.


Reference: `vasp.bands` outputs
-------------------------------

The following output nodes are created upon successful completion:

* ``bands``, type: :py:class:`aiida.orm.nodes.data.array.bands.BandsData`. The calculated band structure of the material.
* ``misc``, type: :py:class:`aiida.orm.nodes.data.dict.Dict`. Dictionary containing the output parameters containing smaller quantities that do not depend on system size.

Depending on the passed inputs to the workchain several outputs might be exposed according to what was defined in the :ref:`vasp_workchain_outputs`.

.. _VASP: https://www.vasp.at
.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/latest/
.. _seekpath: https://github.com/giovannipizzi/seekpath
