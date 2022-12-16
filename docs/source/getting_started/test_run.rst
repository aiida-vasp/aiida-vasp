.. _test_run:

==================
5. Test a VASP run
==================

Finally, following the previous steps we are ready to launch a `VASP`_ calculation of silicon using the standard PBE silicon potential
that you now should have sitting in your database at this point.

#. Please fetch the ``run_vasp_lean.py`` example::

     $ wget https://github.com/aiida-vasp/aiida-vasp/raw/master/examples/run_vasp_lean.py

#. Make sure your `AiiDA`_ virtual environment is activated. Then, check that the
   daemon is running::

     $ verdi daemon status
     Profile: default
     Daemon is running as PID 78693 since 2022-11-23 13:15:54
     Active workers [1]:
     PID    MEM %    CPU %  started
     -----  -------  -------  -------------------
     78696    0.605        0  2022-11-23 13:15:54
     Use verdi daemon [incr | decr] [num] to increase / decrease the amount of workers

   In case it does not run, start it with ``verdi daemon start``.
     
#. Have a look at the ``run_vasp_lean.py`` file and see if you can intuitively understand
   what it does. Most likely, you need to modify the following entries::

     CODE_STRING = 'vasp@mycluster'

   which is a reference to the `code` we added in :ref:`code` and sets which `code` to use for this
   test run. Furthermore::

     POTENTIAL_FAMILY = 'PBE.54'
     POTENTIAL_MAPPING = {'Si': 'Si'}

   where ``POTENTIAL_FAMILY`` is a reference to the name of the potential family we already have uploaded,
   while ``POTENTIAL_MAPPING`` sets the mapping between the symbol in your `StructureData`_. Here, the default
   is specified, but one could for instance replace the mapping with ``{'Si':'Si_GW'}`` to instead use the
   potentials adjusted for GW. Please consult the `VASP potentials`_ documentation. Now we move on to compute
   specific parameters::

     OPTIONS.account = '' # modify it to contain your account, if the cluster needs it
     OPTIONS.qos = '' # modify it to contain your qos, if the cluster needs it
     OPTIONS.queue_name = '' # modify it to contain your queue name, if the cluster needs it
     OPTIONS.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 16} # modify it to relect your cluster
     OPTIONS.max_wallclock_seconds = 3600 # how long do you think you job will run in seconds?
     OPTIONS.max_memory_kb = 10240000 # how much memory do you need per node (not per CPU) in kB?
     
   for example (local systems might deviate) for the PBS, Torque, Slurm and SGE scheduler, you need to modify ``resources`` as follows::

     OPTIONS.resources = {'num_machines': 1, 'tot_num_mpiprocs': 16, 'parallel_env': 'mpi*'}  # for SGE
     OPTIONS.resources = {'num_machines': 1, 'tot_num_mpiprocs': 16}  # for Slurm/PBS/Torque

   where the ``num_machines`` sets the number of nodes and ``tot_num_mpiprocs`` sets the total number of MPI processes you want
   to launch. Please consult the `AiiDA documentation`_ for additional details.

   There is also a ``KMESH`` and an ``INCAR``, but there is no need to change those as we are only doing a test run now.
   Update whatever needs updating and save the file.

#. Submit the `VASP`_ calculation by executing::

     $ python run_vasp_lean.py

#. Check the progress::

     $ verdi process list -a
     PK  Created    Process label    Process State     Process status
     ----  ---------  ---------------  ----------------  -----------------------------------
     718  2h ago     VaspWorkChain    ⏹ Finished [0]
     720  2h ago     VaspCalculation  ⏹ Finished [0]
     
     Total results: 3
     
     Report: last time an entry changed state: 2h ago (at 14:22:10 on 2022-11-23)
     Report: Using 0% of the available daemon worker slots.

   For your case, this might be slightly different.

   `AiiDA`_ relies mainly on the concept of ref:`workchains` which is a
   composition of a setup and teardown of :ref:`calculations` (or calls to
   other :ref:`workchains`).  A workchain can be composed into one or
   multiple `workflows`. A small amount of basic :ref:`workchains` are
   included in `AiiDA-VASP`_. Users are encouraged to develop new, or
   complementing :ref:`workchains` and submitting them to the repository to
   increase the efficiency and functional reach of all `VASP`_ users. It is after all
   a community effort.

   From the process list above, we see that the test calculation launched one :ref:`workchains` and one :ref:`calculations`.
   Both is in state ``Finished`` and without an exit code, which means a successful run. For
   `AiiDA-VASP`_ we would like users to only think about and use the :ref:`workchains` as an entrypoint and not
   the :ref:`calculations`. The latter is there only to facilitate the execution of the `VASP`_
   executable and to prepare the input files and trigger parsing when calculation is done.

#. Let us have a look at the report of the workchain to check in more detail what happened::

     $ verdi process report 718
     2022-11-23 14:19:58 [15 | REPORT]: [718|VaspWorkChain|run_process]: launching VaspCalculation<720> iteration #1
     2022-11-23 14:22:10 [18 | REPORT]: [718|VaspWorkChain|results]: work chain completed after 1 iterations
     2022-11-23 14:22:12 [20 | REPORT]: [718|VaspWorkChain|on_terminated]: cleaned remote folders of calculations: 720

   Here we can see log outputs from the plugin and `AiiDA`_. For this case, we see that the :ref:`vasp_workchain` have launched
   a :ref:`vasp_calculation` and completed. Also, we see that the remote folders, meaning the folders used to run
   your calculation on the cluster is removed since this is was a successful run, i.e. it completed without an exit code.
   In case, an exit code is raised, the remote folders are not cleaned and then one can log into the remote folder
   by issuing in this case ``verdi calcjob gotocomputer 720``. This will not work for this run as it was successful, but
   works if the remote folder is still present.
   
#. Running a calculation without wanting any output is maybe not so interesting, so let us now inspect what is available
   by default. When a :ref:`calculations` or :ref:`workchains` completes, due to the provenance all relevant links to the
   input and output is stored on the node. A node in `AiiDA`_ is typically all datatypes, integral or custom, :ref:`calculations`
   and :ref:`workchains`. Each input and output also a node itself.
   To check links are stored on the :ref:`vasp_workchain`, issue::

     $ verdi process show 718
     Property     Value
     -----------  ------------------------------------
     type         VaspWorkChain
     state        Finished [0]
     pk           718
     uuid         bde72acd-9d4f-4a94-a3c3-d68d0eca16c3
     label
     description
     ctime        2022-11-23 14:19:55.741795+01:00
     mtime        2022-11-23 14:22:10.734169+01:00
     
     Inputs               PK  Type
     -----------------  ----  -------------
     clean_workdir       715  Bool
     code                  4  InstalledCode
     kpoints             707  KpointsData
     max_iterations      714  Int
     options             711  Dict
     parameters          708  Dict
     potential_family    709  Str
     potential_mapping   710  Dict
     settings            712  Dict
     structure           706  StructureData
     verbose             713  Bool
     
     Outputs          PK  Type
     -------------  ----  ----------
     misc            723  Dict
     remote_folder   721  RemoteData
     retrieved       722  FolderData
     
     Caller      PK  Type
     --------  ----  ---------------
     CALL       717  VerifyWorkChain
     
     Called          PK  Type
     ------------  ----  ---------------
     iteration_01   720  VaspCalculation
     
     Log messages
     ---------------------------------------------
     There are 3 log messages for this calculation
     Run 'verdi process report 718' to see them
     

#. Let us now inspect the outputs a bit. In the outputs there for instance a ``misc``.
   This is a container for properties that does not depend on system size. Typically, total energies, band gaps,
   maximum forces, different tensors etc. Let us see what it contains::

     $ verdi data core.dict show 723
     {
	 "maximum_force": 0.0,
	 "maximum_stress": 17.91131059,
	 "notifications": [],
	 "run_stats": {
	     "average_memory_used": null,
	     "elapsed_time": 5.283,
	     "maximum_memory_used": 180792.0,
	     "mem_usage_base": 30000.0,
	     "mem_usage_fftplans": 800.0,
	     "mem_usage_grid": 1216.0,
	     "mem_usage_nonl-proj": 615.0,
	     "mem_usage_one-center": 6.0,
	     "mem_usage_wavefun": 1285.0,
	     "system_time": 0.275,
	     "total_cpu_time_used": 1.236,
	     "user_time": 0.961
	 },
	 "run_status": {
	     "consistent_nelm_breach": false,
	     "contains_nelm_breach": false,
	     "electronic_converged": true,
	     "finished": true,
	     "ionic_converged": null,
	     "last_iteration_index": [
		 1,
		 10
	     ],
	     "nbands": 8,
	     "nelm": 60,
	     "nsw": 0
	 },
	 "total_energies": {
	     "energy_extrapolated": -10.79613809,
	     "energy_extrapolated_electronic": -10.79613809
	 },
	 "version": "6.3.2"
     }

   Here we can see that we get some basic information about the run, most of which should be self explanatory if
   you are a familiar VASP user.

.. note::
   Notice that most of the parsing is disabled by default to make sure new users do not run out of disk space
   due to storing more than they need. You can control the output of the parsing by adjusting the parameters
   related to the parser as defined in :ref:`parsing`.

#. Let us now check how we can access the raw `VASP`_ files and verify what files are kept after a successful run
   by default. For this, it is useful to utilize the ``verdi shell` command, which gives
   you an iPython shell with all `AiiDA`_ related candy preloaded. Load that shell. Then in that shell you load the ``718``
   node and inspect for instance the ``vasp_output`` file, which contains both the ``stdout`` and ``stderr`` like this::

     $ verdi shell
     Python 3.10.8 (main, Nov  1 2022, 14:18:21) [GCC 12.2.0]
     Type 'copyright', 'credits' or 'license' for more information
     IPython 7.34.0 -- An enhanced Interactive Python. Type '?' for help.
     
     In [1]: node = load_node(718)
     
     In [2]: node.outputs.retrieved.list_object_names()
     Out[2]: 
     ['OUTCAR',
     '_scheduler-stderr.txt',
     '_scheduler-stdout.txt',
     'vasp_output',
     'vasprun.xml']
     
     In [3]: node.outputs.retrieved.get_object_content('vasp_output')
     Out[3]: ' running on    1 total cores\n distrk:  each k-point on    1 cores,    1 groups\n distr:  one band on    1 cores,    1 groups\n vasp.6.3.2 27Jun22 (build Nov 14 2022 16:53:56) complex                        \n  \n POSCAR found type information on POSCAR Si\n POSCAR found :  1 types and       2 ions\n Reading from existing POTCAR\n scaLAPACK will be used\n Reading from existing POTCAR\n LDA part: xc-table for Pade appr. of Perdew\n POSCAR, INCAR and KPOINTS ok, starting setup\n FFT: planning ... GRIDC\n FFT: planning ... GRID_SOFT\n FFT: planning ... GRID\n WAVECAR not read\n entering main loop\n       N       E                     dE             d eps       ncg     rms          rms(c)\nDAV:   1    -0.318722213841E+01   -0.31872E+01   -0.17933E+03   576   0.362E+02\nDAV:   2    -0.109114195761E+02   -0.77242E+01   -0.75893E+01   896   0.514E+01\nDAV:   3    -0.109752291049E+02   -0.63810E-01   -0.63810E-01   672   0.564E+00\nDAV:   4    -0.109754105753E+02   -0.18147E-03   -0.18147E-03   880   0.313E-01\nDAV:   5    -0.109754107113E+02   -0.13603E-06   -0.13638E-06   688   0.577E-03    0.612E+00\nDAV:   6    -0.108513820736E+02    0.12403E+00   -0.80386E-02   688   0.147E+00    0.376E+00\nDAV:   7    -0.107949973378E+02    0.56385E-01   -0.15442E-01   704   0.217E+00    0.192E-01\nDAV:   8    -0.107958937809E+02   -0.89644E-03   -0.34625E-03   616   0.419E-01    0.708E-02\nDAV:   9    -0.107961129176E+02   -0.21914E-03   -0.21477E-04   792   0.111E-01    0.607E-02\nDAV:  10    -0.107961380924E+02   -0.25175E-04   -0.27776E-05   368   0.397E-02\n   1 F= -.10796138E+02 E0= -.10796138E+02  d E =0.000000E+00\n writing wavefunctions\n'

   As you might observe if you are familiar with Python, we do not open a file here when accessing ``vasp_output``, but use a method
   to get the content of the `object`. This is a method that gets the content in the backend using whatever method is relevant.
   The abstraction was done to first utilize the internal object storage in `AiiDA`_, but this also gave other benefits, for instance
   to utilize object store solutions with e.g. S3 interfaces etc.
     
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _AiiDA: https://www.aiida.net
.. _AiiDA documentation: https://aiida.readthedocs.io/projects/aiida-core/en/latest/index.html
.. _VASP: https://www.vasp.at
.. _VASP potentials: https://www.vasp.at/wiki/index.php/Available_PAW_potentials
.. _StructureData: https://aiida.readthedocs.io/projects/aiida-core/en/latest/topics/data_types.html
