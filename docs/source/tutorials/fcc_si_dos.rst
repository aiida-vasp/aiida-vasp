.. _tutorial_fcc_si_dos:

=========================
5. Parsing and quantities
=========================

As you have realized, the plugin does not parse and store many quantities by default. This is done
intentionally in order not to fill up disk space without there being a need for it.
However, for most users, you would most likely want to have access to more useful quantities
after a `VASP`_ calculation are finished and this tutorial touches on this to give you
a short introduction on how this is done and where to find more information about this.
Let us continue with the `VASP`_ tutorials. In particular, let us now extract
the `density of states for FCC Si`_ as done in the second tutorial. We will from here
and in the rest of the tutorials not run `VASP`_ in the regular manner and leave those
exercise to the reader. Please do so to get a feel for how `AiiDA-VASP`_ makes
these operations more standardized and reproducible,
in particular when you need to perform many or repeat calculations.

Again, we assume you have completed the :ref:`previous tutorial<tutorial_fcc_si_step4>`.

As you might have noticed we have now populated the database with the silicon structure
multiple times. Let us try instead to load some of the structures that already are present
to save coding, reuse previous results and save data storage.

In :ref:`previous tutorial<tutorial_fcc_si_step3>`, when developing the call script
``run_fcc_si_workchain``,
we set the ``structure.label`` in ``get_structure`` method.
Please inspect this now and notice that for a lattice constant of 3.9, the ``label`` is
set to ``silicon_at_3_9``. Remember also that we cannot modify structures already present in
the database. This goes for all nodes, structure, calculation, workchain or any other node.
This is due to the strict compliance of data provenance. If we want to modify a node, for instance,
here the case of the ``StructureData`` node, we need to create a new
:py:class:`Structure Data<aiida.orm.StructureData>` using the stored entry as an initializer.
In this case we will not modify the structure, so we will simply use the unmodified one.

#. Let us first create a simple call script that uses the :ref:`vasp_workchain`:

   .. literalinclude:: ../../../tutorials/run_fcc_si_dos.py

   Or you can download it::

     $ wget https://github.com/aiida-vasp/aiida-vasp/raw/master/tutorials/run_fcc_si_dos.py

   Where we have modified the input parameters according to the tutorial
   `density of states for FCC Si`_.

   .. note::
      Importantly, we have told the parser to give us the density of states by adding
      ``add_dos: True`` dictionary entry to the ``parser_settings`` key in the ``settings``
      dictionary. In fact, this is the general way to control what quantities to parse, to enable with
      ``add_<quantity>: True`` or ``add_<quantity>: False``, where ``<quantity>`` is a given
      quantity. Please consult :ref:`parsing` for the supported quantities and additional details.

#. Let us now execute ``run_fcc_si_dos.py`` by issuing the command::

     $ python run_fcc_si_dos.py

#. After a while we check the status::

     $ verdi process list -a
       PK  Created    Process label         Process State     Process status
     ----  ---------  --------------------  ----------------  --------------------------------------------------------------------------------------
     2431  6m ago     VaspWorkChain         ⏹ Finished [0]
     2433  6m ago     VaspCalculation       ⏹ Finished [0]

#. Check the content of the ``VaspWorkChain``::

     $ verdi process show 2431
     Property     Value
     -----------  ------------------------------------
     type         VaspWorkChain
     state        Finished [0]
     pk           2431
     uuid         40ce7bd6-cd38-405e-951e-c56251a0cf1b
     label
     description
     ctime        2022-12-22 11:17:25.967623+01:00
     mtime        2022-12-22 11:19:39.816274+01:00

     Inputs               PK  Type
     -----------------  ----  -------------
     clean_workdir      2430  Bool
     code                818  InstalledCode
     kpoints            2422  KpointsData
     max_iterations     2429  Int
     options            2426  Dict
     parameters         2423  Dict
     potential_family   2424  Str
     potential_mapping  2425  Dict
     settings           2427  Dict
     structure          1529  StructureData
     verbose            2428  Bool

     Outputs          PK  Type
     -------------  ----  ----------
     dos            2436  ArrayData
     misc           2437  Dict
     remote_folder  2434  RemoteData
     retrieved      2435  FolderData

     Called          PK  Type
     ------------  ----  ---------------
     iteration_01  2433  VaspCalculation

     Log messages
     ---------------------------------------------
     There are 3 log messages for this calculation
     Run 'verdi process report 2431' to see them

   And as you can see, ``dos`` is listed in the output, so let us quickly inspect it::

     $ verdi data core.array show 2436
     ...

   where we have truncated the output. You can verify both the structure of the stored
   density of states data and its values.

   Now, as you may already know, running with such a dense k-point grid for the initial
   calculation is usually not a good idea. It is more efficient to pre-converge the
   electronic states using a more sparse k-point grid and then restart the calculation
   using a more dense k-point grid when calculating the density of states. Of course
   this can form the foundations of a workflow dedicated to pre-converging results, not
   just for density of states calculations. We will leave this exercise to the user.

   With this we conclude the initial tutorials that follow the `VASP`_ tutorials closely
   and will now continue with a few tutorials related to how to interact with the plugin
   and finally conclude with a few rather specific tutorials that you might find interesting or
   inspiring.

.. _AiiDA: https://www.aiida.net
.. _density of states for FCC Si: https://www.vasp.at/wiki/index.php/Fcc_Si_DOS
.. _VASP: https://www.vasp.at
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
