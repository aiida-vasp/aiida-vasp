############################
Single Calculation Workflows
############################

***********
Description
***********

There are four workflows in this category. They are mainly intended to be used as building blocks in more complex
workflows, however they can also be used as convenience wrappers around their corresponding Calculation plugins.

The workflows are self-documenting in terms of input parameters, Example::

   ScfWorkflow = WorkflowFactory('vasp.scf')
   parameters = ScfWorkflow.get_params_template()
   print parameters

will print out a dictionary of possible input parameters with explaining comments.
It is also possible to write and load parameters to and from JSON files::

   ScfWorkflow.get_template(path='scf_inputs.json')

Will write the same information to a file, which can be edited with the desired values.
It can then be loaded and fed to a workflow::

   import json
   with open('scf_inputs.json) as inputs:
      parameters = json.load(inputs)
   wf = ScfWorkflow(params=parameters)
   wf.start()

**********
Parameters
**********

* :py:class:`ScfWorkflow <aiida.workflows.vasp.scf.ScfWorkflow>`, runs a :doc:`ScfCalculation <../calcs/scf>`
* :py:class:`NscfWorkflow <aiida.workflows.vasp.nscf.NscfWorkflow>`, runs a :doc:`NscfCalculation <../calcs/nscf>`
* :py:class:`ProjectionsWorkflow <aiida.workflows.vasp.projections.ProjectionsWorkflow>`, runs an :doc:`AmnCalculation <../calcs/amn>`
* :py:class:`WannieWannierrWorkflow <aiida.workflows.vasp.wannier.WannierWorkflow>`, runs an :doc:`WannierCalculation <../calcs/wannier>`

*********
Reference
*********

.. automodule:: aiida.workflows.vasp.scf

   .. autoclass:: ScfWorkflow
      :members: get_params_template, get_template

      .. automethod:: start

         Run a ScfCalculation with input nodes constructed
         from the workflow's parameters

      .. automethod:: end

         Checks wether the run was successful in producing
         wavefunctions and chargedensities nodes. Adds the scf calculation
         node to the results if successful.

.. automodule:: aiida.workflows.vasp.nscf

   .. autoclass:: NscfWorkflow
      :members: get_params_template, get_template

      .. automethod:: start

         Runs an NscfCalculation using the given uuid of a
         ScfCalculation.

      .. automethod:: end

         Checks for success of the NscfCalculation.
         If all the expected outputs are present, adds the calculation
         node to the results.

.. automodule:: aiida.workflows.vasp.projections
   
   .. autoclass:: ProjectionsWorkflow
      :members: get_params_template, get_template

      .. automethod:: start

         Runs an AmnCalculation using the given uuid of a 
         NscfCalculation. Modifies wannier_settings from the nscf calculation
         it continues from transparently using an InlineCalculation.

      .. automethod:: end

         Checks for success of the AmnCalculation.
         If all the expected outputs are present, adds the calculation
         node to the results.

.. automodule:: aiida.workflows.vasp.wannier
   .. autoclass:: WannierWorkflow
      :members: get_params_template, get_template

      .. automethod:: start

         Runs a WannierCalculation using the given uuid of a 
         AmnCalculation. Modifies wannier_settings from the amn calculation
         it continues from transparently using an InlineCalculation, if necessary.

      .. automethod:: end

         Checks for success of the WannierCalculation.
         If all the expected outputs are present, adds the calculation
         node to the results.
