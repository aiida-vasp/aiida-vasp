###############
WindowsWorkflow
###############

***********
Description
***********

Runs all the steps required to get a bandstructure and hopping terms from wannier90 as well as an
ab inition Bandstructure from VASP that can be used to compare.
The window parameters for vasp90 are given manually at the start. Several sets of parameters can be given, for each a wannier workflow is run.

::

   User Input (given at the start)
      |
      |-------+------------------------------+-------------+
      |       |                              |             |
      v       v                              v             v
   Scf ---> Nscf -> (wannier_parameters) ---> Projections ---> Wannier ---> Nscf
                                                   |----> Wannier       |
                                                             |          +---> Bands, DOS
                                                             v
                                                    Bands, wannier90_hr.dat

minimal example::

   WindowsWorkflow = WorkflowFactory('vasp.windows')

   num_nodes = 2
   num_ppn = 6
   nbands = 24
   assert(nbands % (num_nodes * num_ppn) == 0) # avoid getting more bands than asking for

   parameters = {
      "kpoints": {
         "mesh": [4, 4, 4],
         "path": [
            ["X", [0., .5, .5], "G", [0., 0., 0.]]
         ]
      },
      "paw_family": "<name>",
      "paw_map": {"In": "In_d", "As": "As"},
      "resources": {
         "num_machines": num_nodes,
         "num_mpiprocs_per_machine": num_ppn
      },
      "wannier_resources": {
         "num_machines": 1,
         "num_mpiprocs_per_machine": 1
      }
      "queue": "<computing queue name>",
      "vasp_code": "<vasp-code@computer>",
      "wannier_code": "<wannier-code@computer>",
      "structure": "<path to cif or poscar file for InAs>",
      "parameters": {
         "nbands": nbands,
         "ediff": 1e-5,
         "gga": "PE",
         "gga_compat": False,
      }
      "wannier_parameters": {
         "dis_num_iter": 100,
         "num_iter": 100
      }
      "projections": [
         "In : s; px; py; pz",
         "As : s; px; py; py"
      ]
      "windows": [
         {"inner": [<min>, <max>], "outer": [<min>, <max>]},
         ...
      ]
   }

   wf = WindowsWorkflow(params=parameters)
   wf.start()

The k-point path for wannier interpolation as well as for the ab-initio band structure
are both set in the key "kpoints". Any kpoint path set in "wannier_parameters" would be overwritten.
Obvious input parameters like LWANNIER90 for VASP or bands_plot for wannier90 are set automatically.

**********
Parameters
**********

Most of the parameters have the same semantics as for the :ref:`single calculation workflows <stepwf-parameters>`

* parameters
* wannier_parameters
* structure
* paw_family
* paw_map
* resources
* queue
* vasp_code
* wannier_code
* description
* extras
* label

The following parameters are specific to this calculation or have changed semantics versus the single calculation workflow parameters.

* kpoints: dict, like for single calculation workflows, except containing a mesh as well as a path specification.
* projections:, list as it would be set in wannier_parameters
* windows: list of dicts, each dict follows the format {'outer': [min, max], 'inner': [min, max]}
* wannier_resources: overrides resources for wannier calculations
* wannier_queue: overrides queue for wannier calculations (it's possible to have VASP and wannier on separate computers)

*********
Reference
*********

.. automodule:: aiida_vasp.workflows.windows

   .. autoclass:: WindowsWorkflow
      :members: get_params_template, get_template, Helper

      .. automethod:: start

         runs an ScfWorkflow, passing on the relevant parameters.

      .. automethod:: get_win

         runs an NscfWorkflow continuing from the previous step's ScfCalculation.
         sets use_wannier automatically.'

      .. automethod:: get_projections

         runs an AmnWorkflow with the results of the previoust step
         to get the wannier projections file.
         sets some wannier_parameters keys automatically.

      .. automethod:: get_tbmodel

         loops over the given windows parameters and starts a WannierWorkflow
         for each set.

      .. automethod:: get_reference_bands

         starts another NscfWorkflow to get a band structure for the same
         path as the wannier interpolations.
         Can be compared kpoint for kpoint to the wannier interpolated bands.

      .. automethod:: make_results

         adds all obtained bands nodes to the workflow's results.
