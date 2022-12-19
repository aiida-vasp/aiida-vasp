.. _tutorial_fcc_si_dos:

===========================
4. FCC Si density of states
===========================

Let us continue with the `VASP`_ tutorials. In particular, let us now extract
the `density of states for FCC Si`_ as done in the second tutorial. We will from here
and in the rest of the tutorials not run `VASP`_ in the regular manner and leave those
exercise to the reader. Please do so, to appreciate how much simpler `AiiDA-VASP`_ makes
these operations, in particular when you need to perform many or repeat calculations.

Again, we assume you have completed the :ref:`previous tutorial<tutorial_fcc_si_step3>`.

As you might have noticed we have now populated the database with the silicon structure
multiple times. Let us try instead to load some of the structures that already are present
to save coding, reuse previous results and save data storage.

In :ref:`previous tutorial<tutorial_fcc_si_step3>`, when developing the call script,
we set the ``label``. Please inspect this. E.g. for a lattice constant of 3.9, the ``label`` is
set to ``silicon_at_3_9``. Remember also that we cannot modify structures already present in
the database. If we want to do this, we need to create a new
:py:class:`Structure Data<aiida.orm.StructureData>` using the stored entry as an initializer.
In this case we will not modify the structure, so we will simply use the unmodified one.

#. Let us first create a simple call script that uses the :ref:`vasp_workchain`:

   .. literalinclude:: ../../../tutorials/run_fcc_si_dos.py

   Or you can download it::

     %/$ wget https://github.com/aiida-vasp/aiida-vasp/raw/develop/tutorials/run_fcc_si_dos.py

   Where we have modified the input parameters according to the tutorial
   `density of states for FCC Si`_. And importantly, we have also told the parser to
   give us the density of states.

#. And execute it::

     %/$ python run_fcc_si_dos.py

#. After a while we check the status::

     %/$ verdi process list -a
     PK  Created    Process label                 Process State      Process status
     ------  ---------  ----------------------------  -----------------  --------------------------------------
     ---------------------
     103820  2m ago     VerifyWorkChain               ⏹ Finished [0]
     103821  2m ago     VaspWorkChain                 ⏹ Finished [0]
     103822  2m ago     VaspCalculation               ⏹ Finished [0]

#. Check the content of the topmost workchain::

     %/$ verdi process show 103820
     (aiida) [efl@efl tutorials]$ verdi process show 103820
     Property       Value
     -------------  ------------------------------------
     type           WorkChainNode
     pk             103820
     uuid           ae0047fe-6351-4585-a60c-ad02dcda8f93
     label
     description
     ctime          2019-10-08 14:30:46.965520+00:00
     mtime          2019-10-08 14:31:53.971080+00:00
     process state  Finished
     exit status    0
     computer       [6] mycluster

     Inputs                     PK  Type
     ---------------------  ------  -------------
     clean_workdir          103818  Bool
     code                   101271  Code
     kpoints                103810  KpointsData
     max_iterations         103817  Int
     options                103814  Dict
     parameters             103811  Dict
     potential_family       103812  Str
     potential_mapping      103813  Dict
     settings               103815  Dict
     structure              103729  StructureData
     verbose                103816  Bool
     verify_max_iterations  103819  Int

     Outputs            PK  Type
     -------------  ------  ----------
     dos            103825  ArrayData
     misc           103826  Dict
     remote_folder  103823  RemoteData
     retrieved      103824  FolderData

     Called        PK  Type
     --------  ------  -------------
     CALL      103821  WorkChainNode

     Log messages
     ---------------------------------------------
     There are 1 log messages for this calculation
     Run 'verdi process report 103820' to see them

   And as you can see, ``dos`` is now listed in the output.

   Now, as you may already know, running with such a dense k-point grid for the initial
   calculation is usually not a good idea. It is more efficient to pre-converge the
   electronic states using a more sparse k-point grid and then restart the calculation
   using a more dense k-point grid when calculating the density of states.

.. _AiiDA: https://www.aiida.net
.. _density of states for FCC Si: https://cms.mpi.univie.ac.at/wiki/index.php/Fcc_Si_DOS
.. _VASP: https://www.vasp.at
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
