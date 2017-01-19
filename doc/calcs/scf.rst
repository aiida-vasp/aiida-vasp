##############
ScfCalculation
##############

***********
Description
***********
This Calculation plugin runs a self-consistent VASP calculation and retrieves the CHGCAR and WAVECAR files needed for running further non-selfconsistent calculations. Another possible application would be to run a series of SCFCalculations to search for the optimal lattice parameter.

******
Inputs
******
* :ref:`settings <vasp-input-settings>`
* :ref:`kpoints <vasp-input-kpoints>`, kpoints must be set via :py:meth:`set_kpoints_mesh <aiida.orm.data.array.kpoints.KpointsData.set_kpoints_mesh>`
* :ref:`structure <vasp-input-structure>`
* :ref:`paw <vasp-input-paw>`, one per element in the structure.

*******
Outputs
*******
* :ref:`results <vasp-output-results>`
* :ref:`kpoints <vasp-output-kpoints>`
* :ref:`chargedensities <vasp-output-chargedens>`
* :ref:`wavefunctions <vasp-output-wavefun>`

*********
Reference
*********
Derived from :py:class:`BasicCalculation <aiida_vasp.calcs.base.BasicCalculation>`

.. autoclass:: aiida_vasp.calcs.scf.ScfCalculation
   :members: verify_inputs, _prepare_for_submission
   :undoc-members:
