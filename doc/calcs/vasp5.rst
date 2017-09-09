################
Vasp5Calculation
################

***********
Description
***********

This calculation retrieves and stores all possible output files it finds after the run, including files obtained using the
Vasp2Wannier90 interface. Therefore it should only be used if every single file is needed for later steps or
if it is unclear which files might be required, in order to avoid storage of unnecessary large files.
Most of the files will be stored in the file repository but not parsed into a database node by default.

Therefore it may be necessary to write an InlineCalculation to take the "retrieved" node as an input and
output data nodes for files that will be used later on. Otherwise the data proveniency chain connecting input parameters,
codes, outpus and analyzed data will be broken.

******
Inputs
******

* :ref:`parameters <vasp-input-parameters>`
* :ref:`kpoints <vasp-input-kpoints>`
* :ref:`structure <vasp-input-structure>`
* :ref:`paw <vasp-input-paw>`
* :ref:`chargedensities <vasp-input-chargedens>`
* :ref:`wavefunctions <vasp-input-wavefunctions>`

*******
Outputs
*******

* :ref:`results <vasp-output-results>`
* :ref:`bands <vasp-output-bands>`
* :ref:`dos <vasp-output-dos>`
* :ref:`wannier_parameters <vasp-output-wannier_parameters>`
* :ref:`wannier_data <vasp-output-wannier_data>`

***************
Retrieved Files
***************

See `Files used by VASP <http://cms.mpi.univie.ac.at/vasp/Files_used_by_VASp.html>` for reference on VASP files, as well as `wannier90 User Guide <http://www.wannier.org/doc/user_guide.pdf>` for explanations about wannier90 files.

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

*********
Reference
*********
Superclasses:

* :py:class:`NscfCalculation <aiida_vasp.calcs.nscf.NscfCalculation>`
* :py:class:`WannierBase <aiida_vasp.calcs.wannier.WannierBase>`

.. autoclass:: aiida_vasp.calcs.vasp5.Vasp5Calculation
   :members: verify_inputs, _prepare_for_submission
   :undoc-members:
