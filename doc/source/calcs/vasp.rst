################
VaspCalculation
################

***********
Description
***********

This calculation type is used to run VASP.

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

***************
Retrieved Files
***************

See `Files used by VASP <http://cms.mpi.univie.ac.at/vasp/Files_used_by_VASp.html>` for reference on VASP files.

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

*********
Reference
*********
Superclasses:

* :py:class:`VaspCalcBase <aiida_vasp.calcs.base.VaspCalcBase>`

.. autoclass:: aiida_vasp.calcs.vasp.VaspCalculation
   :members: verify_inputs, _prepare_for_submission
   :undoc-members:
