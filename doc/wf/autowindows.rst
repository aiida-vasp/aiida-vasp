###################
AutowindowsWorkflow
###################

***********
Description
***********

Very similar to :doc:`WindowsWorkflow <windows>`. While the windows workflow requires manual
setting of inner & outer windows for wannier90, autowindows automatically centers the windows
and varies ther size by given increments.
Outer windows are automatically set to contain minimum num_bands amount of bands.
Innter windows are generated centered around and always including the band gap.
The user then only has to input how many outer and inner window sizes should be used and
by what increment they should be different.

minimal example::

   AutowindowsWorkflow = WorkflowFactory('vasp.autowindows')

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

   wf = AutowindowsWorkflow(params=parameters)
   wf.start() 
 
**********
Parameters
**********

The only difference to :doc:`WindowsWorkflow <windows>` is that the windows parameter is replaced
by parameters controlling the automatic generation of windows.

* num_owindows: int, number of outer windows to generate.
* owindows-increment: float, size difference between any two outer windows.
* num_iwindows: int, number of inner windows to generate.
* iwindows-increment: float, size difference between any two inner windows.

*********
Reference
*********

.. automodule:: aiida.workflows.vasp.autowindows
   
   .. autoclass:: AutowindowsWorkflow
      :members: get_params_template, get_template, Helper

      .. automethod:: start

         runs the ScfWorkflow

      .. automethod:: get_win

         runs a Vasp2Wannier90 type calculation to get wannier input

      .. automethod:: get_projections

         runs a ProjectionsWorkflow to get the wannier projections file

      .. automethod:: get_bands_preview

         runs a low resolution bandstructure calculation NscfWorkflow
         to get the necessary values for choosing wannier90 windows in the next
         step

      .. automethod:: make_windows

         reads window creation parameters, generates window parameters
         and stores them in an attribute of the workflow

      .. automethod:: get_tbmodel

         runs wannier for interpolated bands with the window parameters
         generated in the previous step

      .. automethod:: get_reference_bands

         runs a bandstructure NscfWorkflow for the same kpoints as the
         wannier interpolated bandstructures.

      .. automethod:: gather_results

         adds links to the calculated bandstructures to the workflow's results.
         
