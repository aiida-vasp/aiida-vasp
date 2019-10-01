.. _code:

================
9. Create a code
================

Now we need to add the code (in this case `VASP`_) to `AiiDA`_.  The
``verdi`` subcommand ``code`` makes it possible to set up a code that later
can be used to execute a given calculation. The code has to be installed on
the `computer`, i.e. in this case ``mycluster``. Before configuring the code
in `AiiDA`_, please make sure it runs and functions as normal on ``mycluster``.

#. Let us now add the code, which we will call ``vasp`` as follows::
     
     % verdi code setup
     Info: enter "?" for help
     Label: vasp
     Description []:
     Default calculation input plugin: ?
     Info: Default calculation plugin to use for this code.
     Select one of:
           arithmetic.add
           templatereplacer
           vasp.vasp
           vasp.vasp2w90
     Default calculation input plugin: vasp.vasp
     Installed on target computer? [True]:
     Computer: mycluster
     Remote absolute path: /usr/local/calc/vasp/vasp
     Success: Code<1> vasp@mycluster created

   During the setup, you are asked to enter the prepend and append text.
   This is any command you would run before your code is executed. For most cluster
   systems you would need to load the correct modules in the prepend text. In the append
   text section you typically enter cleanup routines etc.
   
#. We can now check if the code is present by issuing::

   % verdi code list
   # List of configured codes:
   # (use 'verdi code show CODEID' to see the details)
   * pk 1 - vasp@mycluster
   # No codes found matching the specified criteria.

#. And look at its details::

   % verdi code show vasp@mycluster
   --------------------  ------------------------------------
   PK                    1
   UUID                  bafec878-3ca5-4f30-9bb1-0144fb760fa0
   Label                 vasp
   Description
   Default plugin        vasp.vasp
   Type                  remote
   Remote machine        boston
   Remote absolute path  /usr/local/calc/vasp/vasp
   Prepend text          No prepend text
   Append text           No append text
   --------------------  ------------------------------------

.. _VASP: https://www.vasp.at
.. _AiiDA: https://www.aiida.net
