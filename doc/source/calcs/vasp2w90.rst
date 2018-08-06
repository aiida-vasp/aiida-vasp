###################
Vasp2w90Calculation
###################

***********
Description
***********

Works exactly the same as :doc:`VaspCalculation <vasp>`, from which it is derived, with the difference
that it takes a ``wannier_parameters`` input as well, for use with the VASP2Wannier90 interface. This requires a VASP code compiled with ``LWANNIER90=True``.

******
Inputs
******

* :ref:`parameters <vasp-input-parameters>`
* :ref:`kpoints <vasp-input-kpoints>`
* :ref:`structure <vasp-input-structure>`
* :ref:`paw <vasp-input-paw>`
* :ref:`chargedensities <vasp-input-chargedens>`
* :ref:`wavefunctions <vasp-input-wavefunctions>`
* :ref:`wannier_parameters <vasp-input-wannier_parameters>`

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

