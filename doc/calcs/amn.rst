##############
AmnCalculation
##############

***********
Description
***********
Named after the extension (.amn) of the file wannier90 uses to store matrix elements for projections. It's main application is to obtain said file.

Intended for running VASP with the LWANNIER90 input parameter set to True and a wannier_settings node containing a "projections" block. The result is a wannier_data output node containing all the files necessary for running wannier.x with "hr_plot=True".

The user is responsible for making sure that the inputs for VASP are consistent with what's given in wannier_settings. This is guaranteed to be true if the same basic inputs are used for the whole chain of runs leading to the chargedensities, wavefunctions and wannier_settings nodes, while also making sure that the total number of cores requested for each was an integer divisor of NBANDS.

******
Inputs
******
* :ref:`settings <vasp-input-settings>`
* :ref:`kpoints <vasp-input-kpoints>`
* :ref:`structure <vasp-input-structure>`
* :ref:`paw <vasp-input-paw>`, one per element in the structure.
* :ref:`chargedensities <vasp-input-chargedens>`
* :ref:`wavefunctions <vasp-input-wavefunctions>`
* :ref:`wannier_settings <vasp-input-wannier_settings>`

*******
Outputs
*******
* :ref:`results <vasp-output-results>`
* :ref:`wannier_data <vasp-output-wannier_data>`

*********
Reference
*********
Superclasses:

* :py:class:`NscfCalculation <aiida.orm.calculation.job.vasp.nscf.NscfCalculation>`
* :py:class:`WannierBase <aiida.orm.calculation.job.vasp.wannier.WannierBase>`

.. autoclass:: aiida.orm.calculation.job.vasp.amn.AmnCalculation
   :members: verify_inputs, _prepare_for_submission
   :undoc-members:
