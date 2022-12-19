.. _code:

================
3. Create a code
================

Now we need to add the `code` (in this case representing the `VASP`_ executable) to the current `AiiDA`_ profile.  The
``verdi`` subcommand ``code`` makes it possible to set up a code that later
can be used to execute a given calculation. Please consult the `AiiDA documentation`_ for more details. The `VASP`_ code should be installed on
a `computer`, preferably on a HPC cluster. In this example, we assume you have already added a computer with a label ``mycluster`` according to
the `AiiDA documentation`_. Before configuring the code in `AiiDA`_, please make sure it runs and functions as normal on the computer ``mycluster``.

#. Let us now add the code, which we will call ``vasp`` as follows::

     % verdi code create core.code.installed
     Report: enter ? for help.
     Report: enter ! to ignore the default and set no value.
     Computer: mycluster
     Absolute filepath executable: /cluster/software/vasp/vasp6.3.2/vasp
     Label: vasp
     Description: VASP 6.3.2 standard version
     Default `CalcJob` plugin: vasp.vasp
     Escape using double quotes [y/N]:
     Success: Created InstalledCode<6>

   The `Absolute filepath executable` is the full path to the `VASP`_ executable installed on the remote `computer`, here
   labeled ``mycluster``. In many cases, you might need to utilize different versions of a `VASP`_ executable, for instance
   in the gamma only or non-collinear configuration. Or with additional auxiliary libraries, like BEEF included. You need to
   add a dedicated `code` for each `VASP`_ version you want to utilize. Also, note that your identifier, or `PK`, which is here ``6``
   will likely be different for you, depending on how many existing things are present in your `AiiDA`_ profile.

   During the end of the setup, you are asked to enter the prepend and append text.
   The prepend section is any command you would run before your code is executed. For most cluster
   systems you would need to load the correct modules in the prepend text. That is the first
   file we edit. Enter something along the lines of::

     module purge
     module load <myvaspmodule>

   in the first open section of the file. Read the comments if you are unsure where to add it. Here, ``myvaspmodule``
   is the name of the `VASP`_ module you will associate with the addded `code`. The `Absolute filepath executable`
   can be obtain with for instance ``module show <myvaspmodule>`` under the ``PATH`` environment variable. You will have to check
   what the name of the actual `VASP`_ executable is as this could be tailored by your HPC maintainers. The default `VASP`_ build system
   yields the ``vasp_std``, ``vasp_ncl`` and ``vasp_gam``, which is maybe a good start. Save and close the file.
   A new file opens, which is related to the append, e.g. what is done after the executable has been executed. In the append
   text section you typically enter cleanup routines etc. If you have none, just save and close the file.

#. We can now check if the code is present by issuing::

     % verdi code list
     Full label                                                                                     Pk  Entry point
     -------------------------------------------------------------------------------------------  ----  -------------------
     vasp@mycluster                                                                                  6  core.code.installed

     Use `verdi code show IDENTIFIER` to see details for a code

   This might of course be different for you if you have multiple entries or have named things differently than in this
   example.

#. Also, let us inspect look its details::

     % verdi code show 6
     --------------  -------------------------------------
     PK              6
     UUID            e5de4c5a-ca44-4a3b-a42c-e5f7e1c21cbb
     Label           vasp
     Description     VASP 6.3.2 standard version
     Default plugin  vasp.vasp
     Prepend text    module purge
                     module load <myvaspmodule>
     Append text
     --------------  ------------------------------------

.. _VASP: https://www.vasp.at
.. _AiiDA: https://www.aiida.net
.. _AiiDA documentation: https://aiida.readthedocs.io/projects/aiida-core/en/latest/index.html
