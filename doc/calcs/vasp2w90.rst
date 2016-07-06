###################
Vasp2w90Calculation
###################

***********
Description
***********

Works exactly the same as :doc:`Vasp5Calculation <vasp5>`, from which it is derived, with the difference
that it takes a wannier_settings input as well, so it covers every possible use case of VASP with LWANNIER90=True.

******
Inputs
******

* :ref:`settings <vasp-input-settings>`
* :ref:`kpoints <vasp-input-kpoints>`
* :ref:`structure <vasp-input-structure>`
* :ref:`paw <vasp-input-paw>`
* :ref:`chargedensities <vasp-input-chargedens>`
* :ref:`wavefunctions <vasp-input-wavefunctions>`
* :ref:`wannier_settings <vasp-input-wannier_settings>`

*******
Outputs
*******

* :ref:`results <vasp-output-results>`
* :ref:`bands <vasp-output-bands>`
* :ref:`dos <vasp-output-dos>`
* :ref:`wannier_settings <vasp-output-wannier_settings>`
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

* :py:class:`Vasp5Calculation <aiida.orm.calculation.job.vasp.vasp5.Vasp5Calculation>`

.. autoclass:: aiida.orm.calculation.job.vasp.vasp2w90.Vasp2w90Calculation
   :members: verify_inputs, _prepare_for_submission
   :undoc-members:
