Vasp5Calculation
----------------

Applications
++++++++++++

This calculation retrieves and stores all possible output files it finds after the run, including files obtained using the
Vasp2Wannier90 interface. Therefore it should only be used if every single file is needed for later steps or
if it is unclear which files might be required, in order to avoid storage of unnecessary large files.
Most of the files will be stored in the file repository but not parsed into a database node by default.

Therefore it may be necessary to write an InlineCalculation to take the "retrieved" node as an input and
output data nodes for files that will be used later on. Otherwise the data proveniency chain connecting input parameters,
codes, outpus and analyzed data will be broken.

I/O
+++

Inputs:

* settings: :py:class:`ParameterData <aiida.orm.data.structure.StructureData>`, see :ref:`vasp-settings`.
* kpoints: :py:class:`KpointsData <aiida.orm.data.array.kpoints.KpointsData>`, see :ref:`vasp-kpoints`.
* structure: :py:class:`StructureData <aiida.orm.data.structure.KpointsData>` or :py:class:`CifData <aiida.orm.data.cif.CifData>`, 
see :ref:`vasp-structure`, see :ref:`vasp-structure`.
* paw: :py:class:`PawData <aiida.orm.data.vasp.PawData>`, see :ref:`vasp-paw`.

Default Parser: :py:class:`Vasp5Parser <aiida.parsers.plugins.vasp.vasp5.Vasp5Parser>`.

Retrieved Files:

* CHG
* CHGCAR
* CONTCAR
* DOSCAR
* EIGENVAL
* ELFCAR
* IBZKPT
* LOCPOT
* OSZICAR
* OUTCAR
* PCDAT
* PROCAR
* PROOUT
* TMPCAR
* WAVECAR
* XDATACAR
* vasprun.xml
* wannier90* (all files starting with "wannier90")
