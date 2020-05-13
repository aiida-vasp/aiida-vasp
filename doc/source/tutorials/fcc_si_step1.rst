.. _tutorial_fcc_si_step1:

===============================
1. First introduction to FCC Si
===============================

We will now follow the `FCC Si`_ example in the `VASP`_ tutorial. We will
obtain the same data, but using different strategies so that you as a user will
be familiar with what is going on. The `FCC Si`_ tutorial is all about calculating
the total energies at different volumes to find the volume that gives you the lowest
energy.

In this tutorial we will complete the original `FCC Si`_ tutorial without using
`AiiDA-VASP`_ and then continue to execute the first static (fixed volume) calculation
with `AiiDA-VASP`_ in order to get to know how the submission, monitoring and inspection
process works. The :ref:`next tutorial<tutorial_fcc_si_step2>` will start to perform
calculations at different volumes.

#. First complete the `FCC Si`_ tutorial.

#. Enable your `AiiDA`_ virtual environment where `AiiDA-VASP`_ is installed.

#. Make sure your `AiiDA`_ daemon runs::

     %/$ verdi daemon start

#. When the standard tutorial is completed we will now do the same in `AiiDA-VASP`_ using
   different strategies to you get a feel for how you can structurize a simple workflow
   like this. Let us now fetch the `AiiDA-VASP`_ run file for this example::

     %/$ wget https://github.com/aiida-vasp/aiida-vasp/raw/develop/tutorials/run_fcc_si_one_volume.py

#. Inspect the file, which has the following content::

.. literalinclude:: ../../../tutorials/run_fcc_si_one_volume.py

#. Change the following based on what is required by your previous setup steps::

     code_string

   which is the entry you would like to use when executing ``verdi code list``.

#. Change the following to comply with requirements of your cluster or your project::

     options.account = ''
     options.qos = ''
     options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 16}
     options.queue_name = ''

   For example, if you use a SGE scheduler, you need to modify `resources` as follows::
     resources = {'num_machines': 1, 'tot_num_mpiprocs': 16, 'parallel_env': 'mpi*'}  # for SGE

#. Save and execute the resulting run script by issuing::

     %/$ python run_fcc_si_one_volume.py

#. Check its progress with::

     %/$ verdi process list
     PK  Created    Process label                 Process State      Process status
     ------  ---------  ----------------------------  -----------------  -----------------------------------------------------------

     101366  25s ago    VerifyWorkChain               ⏵ Waiting          Waiting for child processes: 101367
     101367  24s ago    VaspWorkChain                 ⏵ Waiting          Waiting for child processes: 101368
     101368  24s ago    VaspCalculation               ⏵ Waiting          Monitoring scheduler: job state RUNNING

   Different messages appear, depending on when you run ``verdi process list``. Notice that only
   active processes are shown by this command. In order to also list the processes that has
   completed we use ``verdi process list -a``.

#. After a while, we execut ``verdi process list -a`` and get::

     %/$ verdi process list -a
     PK  Created    Process label                 Process State      Process status
     ------  ---------  ----------------------------  -----------------  -----------------------------------------------------------

     101366  3m ago     VerifyWorkChain               ⏹ Finished [0]
     101367  3m ago     VaspWorkChain                 ⏹ Finished [0]
     101368  3m ago     VaspCalculation               ⏹ Finished [0]

   The processes composing the workflow are now finished.
   Basically, your calculation, or workflow is done. As you can see,
   the run were composed of three :ref:`workchains`.

#. Let us have a look at what happened (we typically inspect the topmost, i.e. the workchain with the
   lowest PK)::

     %/$ verdi process report 101366
     2019-10-02 11:02:37 [5034 | REPORT]: [101366|VerifyWorkChain|run_next_workchain]: launching VaspWorkChain<101367> iteration #1
     2019-10-02 11:02:37 [5035 | REPORT]:   [101367|VaspWorkChain|run_calculation]: launching VaspCalculation<101368> iteration #1
     2019-10-02 11:03:25 [5036 | REPORT]:   [101367|VaspWorkChain|_handle_succesfull]: VaspCalculation<101368> completed successfully
     2019-10-02 11:03:26 [5037 | REPORT]:   [101367|VaspWorkChain|results]: VaspWorkChain<101367> completed after 1 iterations
     2019-10-02 11:03:26 [5038 | REPORT]:   [101367|VaspWorkChain|results]: attaching the node Dict<101371> as 'misc'
     2019-10-02 11:03:26 [5039 | REPORT]:   [101367|VaspWorkChain|results]: attaching the node RemoteData<101369> as 'remote_folder'
     2019-10-02 11:03:26 [5040 | REPORT]:   [101367|VaspWorkChain|results]: attaching the node FolderData<101370> as 'retrieved'
     2019-10-02 11:03:28 [5041 | REPORT]:   [101367|VaspWorkChain|on_terminated]: cleaned remote folders of calculations: 101368

#. Let us have a look at what is stored on the topmost workchain. The topmost workchain typically
   contain the output of the workflow calculation::

     %/$ verdi process show 101366
     Property       Value
     -------------  ------------------------------------
     type           WorkChainNode
     pk             101366
     uuid           c914f898-cf59-4fad-8ea8-a189db9af379
     label
     description
     ctime          2019-10-02 11:02:36.177448+00:00
     mtime          2019-10-02 11:03:27.471046+00:00
     process state  Finished
     exit status    0
     computer       [6] mycluster

     Inputs                     PK  Type
     ---------------------  ------  -------------
     clean_workdir          101364  Bool
     code                   101271  Code
     kpoints                101356  KpointsData
     max_iterations         101363  Int
     options                101360  Dict
     parameters             101357  Dict
     potential_family       101358  Str
     potential_mapping      101359  Dict
     settings               101361  Dict
     structure              101355  StructureData
     verbose                101362  Bool
     verify_max_iterations  101365  Int

     Outputs            PK  Type
     -------------  ------  ----------
     misc           101371  Dict
     remote_folder  101369  RemoteData
     retrieved      101370  FolderData

     Called        PK  Type
     --------  ------  -------------
     CALL      101367  WorkChainNode

     Log messages
     ---------------------------------------------
     There are 1 log messages for this calculation
     Run 'verdi process report 101366' to see them


   Here you can see the inputs and outputs of your workflow (attached on a workchain).

#. Let us inspect the ``misc`` output::

     %/$ verdi data dict show 101371
     {
         "maximum_force": 0.0,
         "maximum_stress": 20.21457093,
         "symmetries": {
             "num_space_group_operations": {
                  "dynamic": [
                      48
		  ],
		  "static": [
		      48
		  ]
	     },
	     "point_group": {
	         "dynamic": [
                     null
		 ],
		 "static": [
                     null
		 ]
	     },
	     "primitive_translations": [
	         1
	     ]
	 },
	 "total_energies": {
             "energy_no_entropy": -4.87588342
	 }
     }

   As you can see, this contains the ``maximum_force``, ``maximum_stress`` and ``total_energies``
   in standard VASP units. The container ``misc`` is used to house quantities that are not
   system size dependent (with size, we also mean grid sizes etc., like the k-point grid).

.. _AiiDA: https://www.aiida.net
.. _FCC Si: https://cms.mpi.univie.ac.at/wiki/index.php/Fcc_Si
.. _VASP: https://www.vasp.at
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
