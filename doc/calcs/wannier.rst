##################
WannierCalculation
##################

***********
Description
***********

Calculation for the wannier.x program from wannier90. Retrieves all output files but parses only the bandstructure and hopping terms into output nodes, if they are output by wannier.x.
To access all other output files, the 'retrieved' output node common to all calculations can be used.

.. _wannier-input-parameters:
.. _wannier_input-data:

******
Inputs
******

* parameters, :py:class:`ParameterData <aiida.orm.data.paramters.ParameterData>`, see :ref:`wannier_parameters<vasp-input-wannier_parameters>` for more details.
* data, :py:class:`ArchiveData <aiida.orm.data.vasp.archive.ArchiveData>`, see :ref:`wannier_data <vasp-input-wannier_data>` for details.

.. _wannier-output-bands:
.. _wannier-output-tb_model:

*******
Outputs
*******

* bands, :py:class:`BandsData <aiida.orm.data.array.bands.BandsData>`, if parameters included bands_plot = True. Band structure interpolated by wannier90.
* tb_model, :py:class:`SinglefileData <aiida.orm.data.singlefile.SinglefileData>`, containing wannier90_hr.dat if parameters included hr_plot = True.

*********
Reference
*********
Superclasses:

* :py:class:`WannierBase <aiida_vasp.calcs.wannier.WannierBase>`

.. autoclass:: aiida_vasp.calcs.wannier.WannierCalculation
   :members: verify_inputs, _prepare_for_submission
   :undoc-members:
