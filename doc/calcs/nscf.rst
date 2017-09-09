###############
NscfCalculation
###############

***********
Description
***********
The NscfCalculation is used to calculate bandstructures, DOS, as well as obtain wannier90 input files.

Technically the NscfCalculation is more flexible than that and will accept any values for ICHARG in the parameters node, not only values > 10.
If a chargedensity node is given the ICHARG value is not set or set to a value other than 1 or 11, no 
CHGCAR file will be used in the calculation, because VASP would not read from it, and a warning is written to the calculation's log.
This behaviour is consistent with running VASP manually, while being slightly more clear about what is actually happening.

******
Inputs
******
* :ref:`parameters <vasp-input-parameters>`
* :ref:`kpoints <vasp-input-kpoints>`, to run with LWANNIER90, kpoints must be set via :py:meth:`set_kpoints_mesh <aiida.orm.data.array.kpoints.KpointsData.set_kpoints_mesh>`. Otherwise all methods of setting kpoints are allowed.
* :ref:`structure <vasp-input-structure>`
* :ref:`paw <vasp-input-paw>`, one per element in the structure.
* :ref:`chargedensities <vasp-input-chargedens>`, not necessary for certain values of ICHARG
* :ref:`wavefunctions <vasp-input-wavefunctions>`

*******
Outputs
*******
* :ref:`results <vasp-output-results>`
* :ref:`bands <vasp-output-bands>`
* :ref:`dos <vasp-output-dos>`
* :ref:`wannier_parameters <vasp-output-wannier_parameters>`
* :ref:`wannier_data <vasp-output-wannier_data>`

*********
Reference
*********
Derived from :py:class:`BasicCalculation <aiida_vasp.calcs.base.BasicCalculation>`

.. autoclass:: aiida_vasp.calcs.nscf.NscfCalculation
   :members: verify_inputs, _prepare_for_submission, write_chgcar, write_wavecar
   :undoc-members:
