.. _vasp_calculation:

================
VASP calculation
================

Inputs
------

* :ref:`parameters <vasp-input-parameters>`
* :ref:`dynamics <vasp-input-dynamics>`
* :ref:`structure <vasp-input-structure>`
* :ref:`potential <vasp-input-potential>`
* :ref:`kpoints <vasp-input-kpoints>`
* :ref:`charge_density <vasp-input-charge>`
* :ref:`wavefunctions <vasp-input-wave>`

.. _vasp-input-parameters:

parameters
^^^^^^^^^^
The input parameters (the content of the INCAR file). This is of the `AiiDA`_ data type :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`, containing keys-value pairs that would be given in an INCAR file when running `VASP`_ manually. The conversion (parsing) from ``parameters`` to a `VASP`_ INCAR file is handled by :py:class:`Incar<parsevasp.incar.Incar>`. Example::

   parameters = Dict(dict={
      'system': 'System Name',
      'nbands': 24,
      'gga': 'PE',
      'gga_compat: False,
      'encut': 280.0}
   )

Key names are case-insensitive, which should be considered when querying the database. Values can be given in their python representation or as strings (strings also may include unit specifications).

.. _vasp-input-dynamics:

dynamics
^^^^^^^^
A dictionary containing information about how the calculation should obey the dynamics of the system. Currently only one flag is respected, that is the ``positions_dof`` which contains a list, one entry per position which contains another list, e.g. ``[True, True, False]`` that describes the flags to set in ``POSCAR`` after each position in order to control the selective dynamics. Example::

  dynamics = Dict(dict={'positions_dof':
    [[True, True, True],
     [False, True, False]
    ]}
  )

.. _vasp-input-kpoints:

kpoints
^^^^^^^
The k-point mapping of the reciprocal space. This is of the `AiiDA`_ data type :py:class:`KpointsData<aiida.orm.nodes.data.array.kpoints.KpointsData>` and can either be given as a mesh, list or path. This information is transformed into either of two kinds of KPOINTS file; a mesh or an explicit list. The transformation of k-point paths to lists of k-points is left to AiiDA to ensure consistency over codes. Mesh files are written as such because `VASP`_ treats them differently than lists and many use cases do not work with lists. The conversion (parsing) from ``kpoints`` to a `VASP`_ KPOINTS file is handled by :py:class:`Kpoints<parsevasp.kpoints.Kpoints>`. Example::

   k_mesh = KpointsData()
   k_mesh.set_kpoints_mesh([4, 4, 4], offset=[0, 0, 0]) # a fairly sparse mesh

This leads to the following KPOINTS::

   Automatic mesh
   0
   Gamma
   4 4 4
   0 0 0

Whereas::

  my_kpoints = [[0.0, 0.0, 0.0],
                [0.1, 0.1, 0.1],
		...
		]
  my_weights = [1., 2., ...]
  assert(len(my_kpoints) == 10)
  assert(len(my_weights) == 10)
  k_list = KpointsData()
  k_list.set_kpoints(my_kpoints)

leads to::

  Explicit list
  10
  Direct
  0.0 0.0 0.0 1.0
  0.1 0.1 0.1 2.0
  ...


To use a k-point path requires knowledge of the structure beforehand::

  structure = CifData.get_or_create('<path-to-cif-file>')
  k_path = KpointsData()
  k_path.set_cell(structure.get_ase().get_cell())
  k_path.set_kpoints_path(value=[('G', 'M'), ('M', ...), ... ])

This leads to::

  Explicit list
  <Number of AiiDA generated kpoints>
  Direct
  0  0  0  1.0
  ...

Look at the class documentation for :py:class:`KpointsData <aiida.orm.nodes.data.array.kpoints.KpointsData>` for more information on how to influence the generation of kpoints from paths. One can also utilize `SeeK-path`_ to create consitent explicit lists of k-points to be used for band structure extractions. This is demonstrated in :ref:`bands_workchain`.

.. _`SeeK-path`: https://github.com/giovannipizzi/seekpath

.. _vasp-input-structure:

structure
^^^^^^^^^
The structure of the atomic layout. This is of the `AiiDA`_ data type :py:class:`StructureData<aiida.orm.nodes.data.structure.StructureData>` or :py:class:`CifData<aiida.orm.nodes.data.cif.CifData>`. The conversion (parsing) to (from) ``structure`` from (to) a `VASP`_ POSCAR file is handled by :py:class:`Poscar<parsevasp.poscar.Poscar>`.

.. _vasp-input-potential:

potential
^^^^^^^^^
A namespace containing the potentials to use for each element (the POTCAR files). This is of the AiiDA-VASP data type :py:class:`PotcarData<aiida_vasp.data.potcar.PotcarData>`. How to upload `VASP`_ POTCAR files can be found at :ref:`potentials`. Once uploaded they can be obtained as follows::

   # input_structure is for instance InAs
   potcar_mapping = {'In': 'In_d', 'As': 'As'}
   potcars = PotcarData.get_potcars_from_structure(structure=input_structure, family_name='PBE.54', mapping=potcar_mapping)

One POTCAR node must be given to the calculations for each element in the system.
The calculations take responsibility for ordering the elements consistently between POSCAR and POTCAR.

.. _vasp-input-charge:

charge_density
^^^^^^^^^^^^^^
The charge density of the electrons (the CHGCAR file). This is of the AiiDA-VASP data type :py:class:`ChargedensityData<aiida_vasp.data.chargedensity.ChargedensityData>` and contains a CHGCAR file from a previous (self-consistent) run. This input is optional.

.. _vasp-input-wave:

wavefunctions
^^^^^^^^^^^^^
The plane wave coefficients (the WAVECAR file). This is of the AiiDA-VASP data type :py:class:`WavefunData<aiida_vasp.data.wavefun.WavefunData>` containing a WAVECAR (or WAVEDER) file from a previous (self-consistent) run. This input is optional.

.. _vasp-input-wannier_parameters:

wannier_parameters
^^^^^^^^^^^^^^^^^^
:py:class:`Dict<aiida.orm.nodes.data.dict.Dict>` containing information that would be given to Wannier90 in a `VASP`_ run with ``LWANNIER90 = TRUE``.

Keyword parameters are mapped to key-value pairs, begin-end blocks are represented as lists with an entry per line.
Numerical and boolean values can be given as python or string representations of the respective type.
An example::

   wannier_parameters = Dict(dict={
      "num_bands": 24,
      "num_wann": 8,
      "projections": [
         ["In: s; px; py; pz"],
         ["As: s; px; py; pz"]
      ]
   })


Outputs
-------

Each `Calculation`_ in `AiiDA`_ has at least the following two output nodes:

* ``retrieved``: An `AiiDA`_ data type :py:class:`FolderData<aiida.orm.nodes.data.folder.FolderData>`, containing information about the folder in the file repository holding the retrieved files after a run of a `Calculation`_ is completed (e.g. a regular `VASP`_ run). Each successfully completed `VASP`_ calculation will retrieve at least vasprun.xml and typically more files.
* ``remote_folder``: An `AiiDA`_ data type :py:class:`RemoteData<aiida.orm.nodes.data.remote.RemoteData>`, containing infomation about the directory on the remote computer where the `Calculation`_ ran.

In addition to input parameters, a number of `VASP`_ specific output nodes may be generated depending on the specific `Calculation`_.

.. _vasp-output-misc:

misc
^^^^
A dictionary container that houses all system size independent properties. It is of an `AiiDA`_ data type
:py:class:`Dict<aiida.orm.nodes.data.dict.Dict>` and contains the keys for the maximum force, stress and total energies.


.. _vasp-output-kpoints:

kpoints
^^^^^^^
:py:class:`KpointsData <aiida.orm.data.array.kpoints.KpointsData>` containing output k-points read from the output file IBZKPT.
This node contains a list of k-points which can be passed to other codes or used to construct input kpoints for a `VASP`_ calculation with hybrid functionals.

Applies to:

* :py:class:`ScfCalculation <aiida_vasp.calcs.scf.ScfCalculation>`
* :py:class:`VaspCalculation <aiida_vasp.calcs.vasp.VaspCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`

.. _vasp-output-chargedens:

chargedensities
^^^^^^^^^^^^^^^
:py:class:`ChargeDensity <aiida.orm.data.vasp.chargedensity.ChargedensityData>` containing the CHGCAR output file.

Applies to:

* :py:class:`ScfCalculation <aiida_vasp.calcs.scf.ScfCalculation>`
* :py:class:`VaspCalculation <aiida_vasp.calcs.vasp.VaspCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`

.. _vasp-output-wavefun:

wavefunctions
^^^^^^^^^^^^^
:py:class:`ChargedensityData <aiida.orm.data.vasp.wavefun.WavefunData>` containing a WAVECAR file from a previous (self-consistent) run.
This input only applies to :py:class:`NscfCalculations <aiida_vasp.calcs.NscfCalculation` and derivates.

Applies to:

* :py:class:`ScfCalculation <aiida_vasp.calcs.scf.ScfCalculation>`
* :py:class:`VaspCalculation <aiida_vasp.calcs.vasp.VaspCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`

.. _vasp-output-bands:

bands
^^^^^
:py:class:`BandsData <aiida.orm.data.array.bands.BandsData>` containing the bands information read from EIGENVAL and/or vasprun.xml.

Applies to:

* :py:class:`NscfCalculations <aiida_vasp.calcs.NscfCalculation`
* :py:class:`VaspCalculation <aiida_vasp.calcs.vasp.VaspCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`

.. _vasp-output-dos:

dos
^^^
:py:class:`ArrayData <aiida.orm.data.array.ArrayData>` containing the DOS information read from DOSCAR and/or vasprun.xml.

Applies to:

* :py:class:`NscfCalculations <aiida_vasp.calcs.NscfCalculation`
* :py:class:`VaspCalculation <aiida_vasp.calcs.vasp.VaspCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`

.. _vasp-output-wannier_parameters:

wannier_parameters
^^^^^^^^^^^^^^^^^^
:py:class:`Dict <aiida.orm.data.parameter.Dict>`
with a representation of the wannier90.win file generated by the VASP2Wannier90 interface, if LWANNIER90=True was given as
an input parameter.

Applies to:

* :py:class:`NscfCalculations <aiida_vasp.calcs.NscfCalculation`
* :py:class:`AmnCalculation <aiida_vasp.calcs.amn.AmnCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`

.. _vasp-output-wannier_data:

wannier_data
^^^^^^^^^^^^
:py:class:`ArchiveData <aiida.orm.data.vasp.archive.ArchiveData>`, holding a compressed tar archive of the wannier_setup output files.

Applies to:

* :py:class:`NscfCalculations <aiida_vasp.calcs.NscfCalculation`
* :py:class:`AmnCalculation <aiida_vasp.calcs.amn.AmnCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`


.. _Calculation: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/calculations.html
.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp
