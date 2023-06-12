.. _tutorial_fcc_si_step1:

===============
1. Introduction
===============

Here we will introduce the plugin. We will now follow the `FCC Si`_ example in the
`VASP tutorials`_ tutorial. We will
obtain the same data, but using different strategies so that you as a user will
be familiar with what is going on with the plugin. If you are a familiar `VASP`_ user what we will
go through should be well known. The `FCC Si`_ tutorial is all about calculating
the total energies at different volumes to find the volume that gives you the lowest
energy.

In this tutorial we will complete the original `FCC Si`_ tutorial without using
`AiiDA-VASP`_ and then continue to execute the first static (fixed volume) calculation
with `AiiDA-VASP`_ in order to get to know how the submission, monitoring and inspection
process works. In the :ref:`next tutorial<tutorial_fcc_si_step2>` will start to perform
calculations at different volumes.

Before starting, we would like you to get a bit familiar with the concepts of `AiiDA`_, so please
have a look at `AiiDA concepts`_ before continuing.

#. First complete the `FCC Si`_ tutorial without using `AiiDA`_ or the `AiiDA-VASP`_ plugin.

#. Enable your `AiiDA`_ virtual environment where `AiiDA-VASP`_ is installed.

#. Make sure your `AiiDA`_ daemon runs::

     $ verdi daemon start

#. We will now following the lines of the `FCC Si`_ tutorial, but using `AiiDA-VASP`_. In the process
   we will be touching different strategies to you get a feel for how you can structurize a simple workflow.
   Let us fetch the `AiiDA-VASP`_ run file for this example::

     $ wget https://github.com/aiida-vasp/aiida-vasp/raw/master/tutorials/run_fcc_si_one_volume.py

#. Inspect the file, which has the following content:

   .. literalinclude:: ../../../tutorials/run_fcc_si_one_volume.py

#. Change the ``CODE_STRING``  based on the code and computer you have stored. This can be
   inspected with ``verdi code list``.

#. Change the following to comply with requirements of your cluster or your project::

     options.account = ''
     options.qos = ''
     options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 8}
     options.queue_name = ''

   For example, if you use a SGE scheduler, you need to modify `resources` as follows::

     resources = {'num_machines': 1, 'tot_num_mpiprocs': 8, 'parallel_env': 'mpi*'}  # for SGE

   Please consult the documentation on `AiiDA job resources`_ and adjust accordingly.

#. Save and execute the resulting run script by issuing::

     $ python run_fcc_si_one_volume.py

#. Check its progress with::

     $ verdi process list
     PK  Created    Process label    Process State     Process status
     ----  ---------  ---------------  ----------------  -----------------------------------
     883  14s ago    VaspWorkChain    ⏵ Waiting         Waiting for child processes: 885
     885  13s ago    VaspCalculation  ⏵ Waiting         Waiting for transport task: upload

     Total results: 2

     Report: last time an entry changed state: 18s ago (at 14:29:24 on 2022-12-19)
     Report: Using 0% of the available daemon worker slots.

   By running ``verdi process list`` we get a list of all active processes. Depending
   on when you run this command, your ``PK`` and specific information shown might be different,
   but the key observation is that we launched a :ref:`vasp_workchain`, which is the main entrypoint
   for launching a simple VASP calculation. This launches a :ref:`vasp_calculation` which is the
   process in `AiiDA`_ which launches the actual VASP calculation. This is presently in the ``Waiting for transport task: upload``
   process status, meaning, that it is currently uploading results to the computer.

   .. note::
      Notice that the ``verdi process list`` only lists the active processes. Hopefully, after a while
      the launched processes will complete without errors and it will not be visible with ``verdi process list``
      any more as finished processes are not considered active. In order to also list these processes we can
      use ``verdi process list -a``.

#. After a while, we execute ``verdi process list -a`` and get::

     $ verdi process list -a
       PK  Created    Process label    Process State     Process status
     ----  ---------  ---------------  ----------------  -----------------------------------
     883  8m ago     VaspWorkChain    ⏹ Finished [0]
     885  8m ago     VaspCalculation  ⏹ Finished [0]

     Total results: 2

     Report: last time an entry changed state: 6m ago (at 17:07:04 on 2022-12-19)
     Report: Using 0% of the available daemon worker slots.

   The processes composing the workflow are now finished. And there is a zero inside the brackets.
   This shows the exit code, and usual practice is that a zero is a sign of a successfully process
   execution. Please consult the documentation of the `AiiDA exit codes`_ for more details. `AiiDA`_
   defines a few internal exit codes and the `AiiDA-VASP`_ plugin adds to those. Consults the documentation
   concerning the specific :ref:`exit_codes`.

   From the finished state we can conclude that your VASP calculation, or workflow is done.

#. Let us have a look at what happened during the execution of the workflow. We typically inspect the topmost,
   i.e. the workchain with the lowest PK as a starting point, here ``883``. Let us look at logs, or report::

     $ verdi process report 883
     2022-12-19 17:04:54 [85 | REPORT]: [883|VaspWorkChain|run_process]: launching VaspCalculation<885> iteration #1
     2022-12-19 17:07:04 [86 | REPORT]: [883|VaspWorkChain|results]: work chain completed after 1 iterations
     2022-12-19 17:07:06 [87 | REPORT]: [883|VaspWorkChain|on_terminated]: cleaned remote folders of calculations: 885

   Nothing particularly interesting and what you would expect.

   .. note::
      Notice that the logs states that the remote folders was cleaned. The default setting of the plugin is to,
      after the :ref:`vasp_workchain` is finished with a zero exit code to clean the remote folder. The remote
      folder is the folder on the computer running the calculations. Typically this is the remote configured cluster
      for VASP calculations. Consult the documentation of :ref:`vasp_workchain` how to modify this behavior if you
      want to change the default setting.

#. Let us have a look at what is stored on the :ref:`vasp_workchain` with ``PK`` of ``883``.
   The topmost workchain typically contain the relevant output of the workflow calculation::

     $ verdi process show 883
     Property     Value
     -----------  ------------------------------------
     type         VaspWorkChain
     state        Finished [0]
     pk           883
     uuid         0c769ee8-07dc-410b-b1eb-7975ca7e7029
     label
     description
     ctime        2022-12-19 17:04:53.027011+01:00
     mtime        2022-12-19 17:07:04.374171+01:00

     Inputs               PK  Type
     -----------------  ----  -------------
     clean_workdir       882  Bool
     code                818  InstalledCode
     kpoints             874  KpointsData
     max_iterations      881  Int
     options             878  Dict
     parameters          875  Dict
     potential_family    876  Str
     potential_mapping   877  Dict
     settings            879  Dict
     structure           873  StructureData
     verbose             880  Bool

     Outputs          PK  Type
     -------------  ----  ----------
     misc            888  Dict
     remote_folder   886  RemoteData
     retrieved       887  FolderData

     Called          PK  Type
     ------------  ----  ---------------
     iteration_01   885  VaspCalculation

     Log messages
     ---------------------------------------------
     There are 3 log messages for this calculation
     Run 'verdi process report 883' to see them

   Here you can see the inputs and outputs of your workflow, which is attached as outputs on a workchain. You can also
   observe the inputs and what other processes have been called, or called this process.

   .. note::
      Most things in `AiiDA`_ that are stored are considered a `node` and we will continue to use this terminology and
      it does not matter if this is an input, output, or a process node, like ``VaspWorkChain``. If you see a ``PK``
      it is for sure a `node`. Please, at this point, reconsider if you need to fresh up on the concepts of
      `AiiDA`_ as explained in `AiiDA concepts`_.

   .. note::
      Most `nodes` can after being stored, which typically is the case when it is passed
      to or from a process, like a workchain or the special calculation process ``VaspCalculation`` not be modified. This
      is a natural consequence of honoring the data provenance concept. At first this can seem a bit frustrating, i.e.
      if you define a `computer`, which is also considered a `node`, do some calculations with this `computer`
      and find out you have to change
      it, you cannot. You have to create a new `computer` with the modified settings. After a while this will come as
      second nature, but takes a bit of getting used to. In the end, if you want data provenance, there is really no
      other good alternative to this.

   .. note::
      Notice that there are three outputs. The ``remote_folder`` gives the path of the remote folder (which is cleaned
      by default if the workflow is considered successful), the ``retrieved``, which is the folder containing the
      retrieved and kept VASP files and ``misc``. For the :ref:`vasp_workchain`, this is the default. If you want to modify
      what is attached in the output, please consult the documentation on :ref:`parsing`.

#. Let us inspect the ``misc`` output::

     $ verdi data core.dict show 888
     {
         "maximum_force": 0.0,
	 "maximum_stress": 20.22402923,
	 "notifications": [],
	 "run_stats": {
	     "average_memory_used": null,
	     "elapsed_time": 2.547,
	     "maximum_memory_used": 47700.0,
	     "mem_usage_base": 30000.0,
	     "mem_usage_fftplans": 296.0,
	     "mem_usage_grid": 451.0,
	     "mem_usage_nonl-proj": 493.0,
	     "mem_usage_one-center": 3.0,
	     "mem_usage_wavefun": 779.0,
	     "system_time": 0.22,
	     "total_cpu_time_used": 0.827,
	     "user_time": 0.607
	 },
	 "run_status": {
	     "consistent_nelm_breach": false,
	     "contains_nelm_breach": false,
             "electronic_converged": true,
             "finished": true,
             "ionic_converged": null,
             "last_iteration_index": [
                 1,
		 8
	     ],
             "nbands": 6,
             "nelm": 60,
             "nsw": 0
	 },
	 "total_energies": {
	     "energy_extrapolated": -4.87588357,
             "energy_extrapolated_electronic": -4.87588357
	 },
	 "version": "6.3.2"
     }

   As you can see, this contains the ``maximum_force``, ``maximum_stress`` and ``total_energies``
   in standard VASP units. The container ``misc`` is used to house quantities that are not
   system size dependent (with size, we also mean grid sizes etc., like the k-point grid, or number of atoms).
   In ``misc`` you also have access to useful run time statistics (mainly what is printed in ``OUTCAR``) and
   run status data.

.. _AiiDA: https://www.aiida.net
.. _FCC Si: https://www.vasp.at/wiki/index.php/Fcc_Si
.. _VASP: https://www.vasp.at
.. _VASP tutorials: https://www.vasp.at/wiki/index.php/Category:Tutorials
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _AiiDA job resources: https://aiida.readthedocs.io/projects/aiida-core/en/latest/topics/schedulers.html?highlight=resources#job-resources
.. _AiiDA exit codes: https://aiida.readthedocs.io/projects/aiida-core/en/latest/topics/processes/usage.html#exit-codes
.. _AiiDA concepts: https://aiida.readthedocs.io/projects/aiida-core/en/latest/topics/provenance/concepts.html
