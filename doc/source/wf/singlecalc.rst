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

For working examples for all workflows, take a look in the vasp plugin repository in
the "examples" folder.

A minimal example::

   ScfWorkflow = WorkflowFactory('vasp.scf')

   num_nodes = 2
   num_ppn = 6
   nbands = 24
   assert(nbands % (num_nodes * num_ppn) == 0) # avoid getting more bands than asking for

   parameters = {
      "kpoints": {"mesh": [4, 4, 4]},
      "paw_family": "<name>",
      "paw_map": {"In": "In_d", "As": "As"},
      "resources": {
         "num_machines": num_nodes,
         "num_mpiprocs_per_machine": num_ppn
      },
      "queue": "<computing queue name>",
      "vasp_code": "<code@computer>",
      "structure": "<path to cif or poscar file for InAs>",
      "parameters": {
         "nbands": nbands,
         "ediff": 1e-5,
         "gga": "PE",
         "gga_compat": false,
      }
   }

   wf = ScfWorkflow(params=parameters)
   wf.start()

.. _stepwf-parameters:

**********
Parameters
**********

* continue_from: string, uuid of an appropriate calculation to continue from.
* parameters: dict, is wrapped into a :py:class:`ParameterData <aiida.orm.data.parameters.ParameterData>` and passed to the calculation. Optional when continuing with a VASP calculation from a previous vasp calculation (can be used to add / orverride keys).
* use_wannier: bool, switches on LWANNIER90 as well as checking for wannier output files in results.
* wannier_parameters: dict, analog to parameters for the wannier_parameters input parameter.
* structure: path to a .cif or POSCAR file. Ignored when continuing from a previous calculation.
* kpoints: dict, containing one and only one of the following keys
   - mesh: list with #kpoints in each direction
   - list: list of kpoints, where each kpoint is a list of coordinates
   - path: list, according to :py:meth:`KpointsData.set_kpoints_path <aiida.orm.data.array.kpoints.KpointsData.set_kpoints_path>`
* paw_family: string, name of a PAW family. Ignored if continuing.
* paw_map: dict, mapping chemical element to PAW symbol. Ignored if continuing.
* queue: string, same as for :py:meth:`JobCalculation.set_queue_name <aiida.orm.calculation.job.JobCalculation.set_queue_name>`, optional when continuing.
* resources: dict, num_machines and num_mpiprocs_per_machine keys as for :py:meth:`JobCalculation.set_resources <aiida.orm.calculation.job.JobCalculation.set_resources>`, optional when continuing.
* vasp_code or wannier_code: string, like for invoking :py:meth:`Code.get_from_string <aiida.orm.code.Code.get_from_string>`.

The following set of parameters can be used to label and categorize the calculations run by the workflow:
* description: string, used to set the description of the calculation.
* extras: dict, passed to :py:meth:`Node.set_extras <aiida.orm.node.Node.set_extras>` of the caclculation.
* label: string, used to set the calculation's label

These properties can be used to filter queries for calculations and therefore to make it easier to find them later in the database.

*********
Workflows
*********

* :py:class:`ScfWorkflow <aiida_vasp.workflows.legacy.scf.ScfWorkflow>`, Obsolete
* :py:class:`NscfWorkflow <aiida_vasp.workflows.legacy.nscf.NscfWorkflow>`, Obsolete
* :py:class:`ProjectionsWorkflow <aiida_vasp.workflows.legacy.projections.ProjectionsWorkflow>`, Obsolete
* :py:class:`WannieWannierrWorkflow <aiida_vasp.workflows.legacy.wannier.WannierWorkflow>`, Obsolete

