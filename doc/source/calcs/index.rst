############
VASP-Plugins
############

***********
Description
***********

There are two types of Calculations for VASP, an all-purpose one for pure VASP calculations and a specialized one for using the VASP2Wannier90 interface.

***************
Supported Codes
***************

Tested with VASP 5.3.3, 5.3.5, 5.4.1.
Versions prior to 5.3.x might behave slightly differently for the same input parameters and output might be formatted differently than
the parsers expect. Some input parameters do not exist in older VASP versions. Therefore interoperability with versions prior to VASP 5.3.x is not guaranteed.

**************
Plugin Details
**************

.. toctree::
   :maxdepth: 4

   All-purpose VASP <vasp>
   Vasp2w90 <vasp2w90>

******
Inputs
******

.. _vasp-input-parameters:

parameters
==========
:py:class:`ParameterData <aiida.orm.data.structure.StructureData>`,
containing keys-value pairs that would be given in an INCAR file when running VASP manually. Example::

   parameters = ParameterData(dict={
      'system': 'System Name',
      'nbands': 24,
      'gga': 'PE',
      'gga_compat: False,
      'encut': 280.0}
   )

Key names are case-insensitive, this should be considered when querying for them.
Values can be given in their python representation or as strings (strings also may include unit specifications).

.. _vasp-input-kpoints:

kpoints
=======
:py:class:`KpointsData <aiida.orm.data.array.kpoints.KpointsData>`,
can either be given as a mesh, list or path. This information is transformed into either of two kinds of KPOINTS file,
mesh or explicit list. The transformation of k-point paths to lists of k-points is left to AiiDA to ensure consistency
over codes. Mesh files are written as such because VASP treats them differently than lists and many use cases do not work with lists.
Examples::

   k_mesh = KpointsData()
   k_mesh.set_kpoints_mesh([4, 4, 4], offset=[0, 0, 0]) # a fairly sparse mesh

This leads to the following KPOINTS::

   Automatic mesh
   0
   Gamma
   4 4 4
   0 0 0

Whereas::

   my_kpoints = [
                 [0, 0, 0],
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
   0   0   0   1.0
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

Look at the class documentation for :py:class:`KpointsData <aiida.orm.data.array.kpoints.KpointsData>` for more information on
how to influence the generation of kpoints from paths.

.. _vasp-input-structure:

structure
=========
:py:class:`StructureData <aiida.orm.data.structure.StructureData>` or :py:class:`CifData <aiida.orm.data.cif.CifData>`.
Use ASE, pymatgen python packages to convert POSCAR files into structure data::

   from ase.io.vasp import read_vasp
   import os
   pwd = os.path.abspath(os.curdir)
   os.chdir(os.path.dirname('<POSCAR absolute path>'))
   atoms = read_vasp('POSCAR')
   os.chdir(pwd)
   structure = self.calc_cls.new_structure()
   structure.set_ase(atoms)

.. _vasp-input-paw:

potcar
======
:py:class:`PawData <aiida.orm.data.vasp.potcar.PotcarData>`, containing POTCAR files.
VASP's POTPAW folders can be uploaded to the database using :: 
   
   verdi data vasp-potcar uploadfamily.

Once uploaded they can be obtained as follows::

   # input_structure is InAs
   potcar_mapping = {'In': 'In_d', 'As': 'As'}
   potcars = PotcarData.load_paw(family='PBE', structure=input_structure, mapping=potcar_mapping)

If for example ``potpaw_PBE/`` was uploaded with the family name "PBE".

One POTCAR input node must be given to the calculations for each element in the system.
The calculations take responsibility for ordering the elements consistently between POSCAR and POTCAR.

.. _vasp-input-chargedens:

charge_density
==============
:py:class:`ChargedensityData <aiida.orm.data.vasp.chargedensity.ChargedensityData>` containing a CHGCAR file from a previous (self-consistent) run.
This input is optional.

.. _vasp-input-wavefunctions:

wavefunctions
==============
:py:class:`ChargedensityData <aiida.orm.data.vasp.wavefun.WavefunData>` containing a WAVECAR file from a previous (self-consistent) run.
This input is optional.

.. _vasp-input-wannier_parameters:

wannier_parameters
==================
:py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>`
containing information that would be given to Wannier90 in a VASP run with LWANNIER90 = TRUE.

Keyword parameters are mapped to key-value pairs, begin-end blocks are represented as lists with an entry per line.
Numerical and boolean values can be given as python or string representations of the respective type.
An example::

   wannier_parameters = ParameterData(dict={
      "num_bands": 24,
      "num_wann": 8,
      "projections": [
         ["In: s; px; py; pz"],
         ["As: s; px; py; pz"]
      ]
   })

.. _vasp-input-wannier_data:

*******
Outputs
*******

Each Calculation in AiiDA has at least the following two output nodes:
* retrieved: :py:class:`FolderData <aiida.orm.data.folder.FolderData>`, containing information about the folder in the file repository holding the retrieved files. Each successfully completed VASP calculation will retrieve at least the OUTCAR, typically more files.
* remote_folder: :py:class:`RemoteData <aiida.orm.data.remote.RemoteData>`, containing info about the folder on the remote computer the calculation ran on.

In addition and depending on the specific Calculation and it's input parameters, a number of VASP-specific output nodes may be generated.

.. _vasp-output-results:

results
=======
:py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>`
containing at least one key 'efermi' with the according value read from vasprun.xml output.
This output node is a good way to store extra results in custom extension calculations.

Applies to all VASP calculations

.. _vasp-output-kpoints:

kpoints
=======
:py:class:`KpointsData <aiida.orm.data.array.kpoints.KpointsData>` containing output k-points read from the output file IBZKPT.
This node contains a list of k-points which can be passed to other codes or used to construct input kpoints for a VASP calculation with hybrid functionals.

Applies to:
* :py:class:`ScfCalculation <aiida_vasp.calcs.scf.ScfCalculation>`
* :py:class:`VaspCalculation <aiida_vasp.calcs.vasp.VaspCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`

.. _vasp-output-chargedens:

chargedensities
===============
py:class:`ChargeDensity <aiida.orm.data.vasp.chargedensity.ChargedensityData>` containing the CHGCAR output file.

Applies to:
* :py:class:`ScfCalculation <aiida_vasp.calcs.scf.ScfCalculation>`
* :py:class:`VaspCalculation <aiida_vasp.calcs.vasp.VaspCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`

.. _vasp-output-wavefun:

wavefunctions
=============
:py:class:`ChargedensityData <aiida.orm.data.vasp.wavefun.WavefunData>` containing a WAVECAR file from a previous (self-consistent) run.
This input only applies to :py:class:`NscfCalculations <aiida_vasp.calcs.NscfCalculation` and derivates.

Applies to:
* :py:class:`ScfCalculation <aiida_vasp.calcs.scf.ScfCalculation>`
* :py:class:`VaspCalculation <aiida_vasp.calcs.vasp.VaspCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`

.. _vasp-output-bands:

bands
=====
:py:class:`BandsData <aiida.orm.data.array.bands.BandsData>` containing the bands information read from EIGENVAL and/or vasprun.xml.

Applies to:
* :py:class:`NscfCalculations <aiida_vasp.calcs.NscfCalculation`
* :py:class:`VaspCalculation <aiida_vasp.calcs.vasp.VaspCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`

.. _vasp-output-dos:

dos
===
:py:class:`ArrayData <aiida.orm.data.array.ArrayData>` containing the DOS information read from DOSCAR and/or vasprun.xml.

Applies to:
* :py:class:`NscfCalculations <aiida_vasp.calcs.NscfCalculation`
* :py:class:`VaspCalculation <aiida_vasp.calcs.vasp.VaspCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`

.. _vasp-output-wannier_parameters:

wannier_parameters
==================
:py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>`
with a representation of the wannier90.win file generated by the VASP2Wannier90 interface, if LWANNIER90=True was given as
an input parameter.

Applies to:
* :py:class:`NscfCalculations <aiida_vasp.calcs.NscfCalculation`
* :py:class:`AmnCalculation <aiida_vasp.calcs.amn.AmnCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`

.. _vasp-output-wannier_data:

wannier_data
============
:py:class:`ArchiveData <aiida.orm.data.vasp.archive.ArchiveData>`, holding a compressed tar archive of the wannier_setup output files.

Applies to:
* :py:class:`NscfCalculations <aiida_vasp.calcs.NscfCalculation`
* :py:class:`AmnCalculation <aiida_vasp.calcs.amn.AmnCalculation>`
* :py:class:`Vasp2w90Calculation <aiida_vasp.calcs.vasp.VaspCalculation>`
