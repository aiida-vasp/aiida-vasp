Kick-start of AiiDA-VASP calculation
=======================================

AiiDA-core and AiiDA-VASP may be installed on your home
directory. Because AiiDA-core and AiiDA-VASP are still under the acitive
development, they are updated frequently and you may be needed to test
across their different versions. So it is a good idea to manage your
AiiDA systems with different versions or settings separately. This
becomes possible to use pyenv or conda environment functionality. Here
the conda environment is used to explain to establish creating your
AiiDA system. The reason to use the later is that it is easy to manage
to install PostgreSQL separately, too.

In the following, the preparation of AiiDA & AiiDA-VASP system is
explained step-by-step. Here, the installation on linux sysytem is
assumed.

Preparing environments
----------------------

Each of AiiDA settings is isolated by first creating its directory. In
it, AiiDA-core, AiiDA-VASP plugin, and PostgreSQL database are
installed as follows. To isolate, ``myaiida`` is used as the
directory, conda environment, and AiiDA setup names.

::

   % cd ~
   % mkdir myaiida
   % cd myaiida
   % conda create -n myaiida python=2.7
   % source activate myaiida
   % conda install gcc_linux-64 gxx_linux-64 ipython

``source activate myaiida`` enables us to enter an independent conda
environment. All the setting or running are always done in this conda
environment. It is mandatory to activate this every time logging to
this computer.

Next setting up PostgreSQL. Below, ``mypassword`` is used as the
PostgreSQL user password.

::

   % conda install postgresql
   % initdb -D /home/username/myaiida
   % pg_ctl -D /home/username/myaiida start
   % createuser -P aiida
   % createdb -O aiida aiidadb
   % psql aiidadb
   aiidadb=# GRANT ALL PRIVILEGES ON DATABASE aiidadb to aiida;
   aiidadb=# \q

The last two lines are executed in the PostgreSQL interactive terminal
invoked by `psql` command.

Install AiiDA and setup profile
-------------------------------

It is ready to install AiiDA after the PostgreSQL setting. Here we use
AiiDA release of `v1.0.0a4`.

::

   % git clone https://github.com/aiidateam/aiida_core
   % cd aiida_core
   % git checkout v1.0.0a4
   % pip install -e .

Here AiiDA setup is done interactively as follows::

   % verdi setup myaiida
   Email Address (identifies your data when sharing): mymail@address.com
   First name: MyFirstName
   Last name: MyLastName
   Institution: Kyoto university
   Setting up profile myaiida.
   AiiDA backend (available: django, sqlalchemy - sqlalchemy is in beta mode): django
   Default user email: mymail@address.com
   Database engine (available: postgresql_psycopg2, mysql - mysql is deprecated): postgresql_psycopg2
   PostgreSQL host: localhost
   PostgreSQL port: 5432
   AiiDA Database name: aiidadb
   AiiDA Database user: aiida
   AiiDA Database password: mypassword
   AiiDA repository directory: /home/username/myaiida/.aiida/repository-myaiida/

To use AiiDA verdi command with tab completion, the file
``$CONDA_PREFIX/etc/conda/activate.d/env_vars.sh`` is created with the
content::

   #!/bin/sh
   export AIIDA_PATH=~/myaiida
   eval "$(_VERDI_COMPLETE=source-bash verdi)"
   export HOST=`hostname`

and ``$CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh``::

   #!/bin/sh
   export HOST=`hostname`

To import this shell environment variables,

::

   % source deactivate
   % source activate myaiida

Now it should be able to start AiiDA daemon by

::

   % verdi daemon start

However it can fail. If failed, run ``verdi setup myaiida`` and
``verdi daemon start`` again.


Install AiiDA-VASP plugin
-------------------------

AiiDA-VASP plugin doesn't have configuration. Installation is done as follows::

   % cd ~/myaiida
   % git clone https://github.com/aiida-vasp/aiida-vasp.git
   % cd aiida-vasp
   % git rev-parse HEAD
   0d7f8569dd4fda298e1ed4307f388428e175e994
   % pip install -e .
   % reentry scan -r aiida


AiiDA's computer setup
----------------------

``computer`` is the computer where our calculations run by the order
of AiiDA daemon. ``computer`` setting includes which computers (e.g.,
cluater, super computer or localhost) are to be used and how our AiiDA
daemon can connect to them.

::

   % verdi computer setup
   Info: enter "?" for help
   Computer label: mycluster
   Hostname: mycluster
   Description []:
   Enable the computer? [True]:
   Transport plugin: ssh
   Scheduler plugin: torque
   Shebang line (first line of each script, starting with #!) [#!/bin/bash]: #!/bin/bash
   Work directory on the computer [/scratch/{username}/aiida/]: /home/username/aiida/
   Mpirun command [mpirun -np {tot_num_mpiprocs}]: mpirun -np {tot_num_mpiprocs}
   Success: Computer<1> mycluster created
   Info: Note: before the computer can be used, it has to be configured with the command:
   Info:   verdi computer configure ssh mycluster

The detail of ssh setting is done by

::

   % verdi computer configure ssh mycluster
   Info: enter "?" for help
   User name [username]:
   port Nr [22]:
   Look for keys [False]:
   SSH key file []: /home/username/.ssh/id_rsa
   Connection timeout in s [60]:
   Allow ssh agent [False]:
   SSH proxy command []:
   Compress file transfers [True]:
   GSS auth [False]:
   GSS kex [False]:
   GSS deleg_creds [False]:
   GSS host [mycluster]:
   Load system host keys [True]:
   Key policy [RejectPolicy]: ?
   Info: SSH key policy
   Select one of:
        RejectPolicy
        WarningPolicy
        AutoAddPolicy
   Key policy [RejectPolicy]: WarningPolicy
   Connection cooldown time (sec) [5]:
   Info: Configuring computer mycluster for user mymail@address.com.
   Success: mycluster successfully configured for mymail@address.com

Test the computer setup by

::

   % verdi computer test mycluster
   Testing computer 'mycluster' for user mymail@address.com...
   > Testing connection...
   > Checking that no spurious output is present...
         [OK]
   > Getting job list...
     `-> OK, 0 jobs found in the queue.
   > Creating a temporary file in the work directory...
     `-> Getting the remote user name...
         [remote username: username]
         [Checking/creating work directory: /home/username/aiida/]
     `-> Creating the file tmpXmpo4J...
     `-> Checking if the file has been created...
         [OK]
     `-> Retrieving the file and checking its content...
         [Retrieved]
         [Content OK]
     `-> Removing the file...
     [Deleted successfully]
   Test completed (all 4 tests succeeded)


AiiDA's code setup
------------------

``code`` describes by which code our calculations run. The code has to
be installed on the location of ``computer``, i.e., if it is a
computer cluster, the code has to be installed properly to run
there. The setup is done as follows::

   % verdi code setup
   Info: enter "?" for help
   Label: vasp544mpi
   Description []:
   Default calculation input plugin: ?
   Info: Default calculation plugin to use for this code.
   Select one of:
        calculation
        function
        inline
        job
        simpleplugins.arithmetic.add
        simpleplugins.templatereplacer
        vasp.vasp
        vasp.vasp2w90
        work
   Default calculation input plugin: vasp.vasp
   Installed on target computer? [True]:
   Computer: mycluster
   Remote absolute path: /usr/local/calc/vasp/vasp544mpi
   Success: Code<1> vasp544mpi@mycluster created

::

   % verdi code list
   # List of configured codes:
   # (use 'verdi code show CODEID' to see the details)
   * pk 1 - vasp544mpi@mycluster
   # No codes found matching the specified criteria.

::

   % verdi code show vasp544mpi@mycluster
   --------------------  ------------------------------------
   PK                    1
   UUID                  bafec878-3ca5-4f30-9bb1-0144fb760fa0
   Label                 vasp544mpi
   Description
   Default plugin        vasp.vasp
   Type                  remote
   Remote machine        boston
   Remote absolute path  /usr/local/calc/vasp/vasp544mpi
   Prepend text          No prepend text
   Append text           No append text
   --------------------  ------------------------------------



Upload PAW dataset to AiiDA database
------------------------------------

To run VASP calculation, PAW datasets have to be written into ``POTCAR``
file. This is generated by AiiDA-VASP plugin. For this, PAW datasets
are stored in AiiDA database. The procedure is as follows::

   % verdi data vasp-potcar uploadfamily --path=/home/username/potpaw_PBE.54.tar --name=PBE.54 --description="PBE potentials for version 5.4"
   skipping file /home/username/potpaw_PBE.54/H_AE/POTCAR - uploading raised <type 'exceptions.IndexError'>list index out of range
   POTCAR files found: 327. New files uploaded: 326, Added to Family: 326

Here it is assumed that the PBE.54 package of the PAW datasets is put
at ``/home/username/potpaw_PBE.54.tar`` as a tar archive.


Run an AiiDA-VASP calculation
-----------------------------

An example of a workchain calculation is copied from ``example``
directory.

   % cd ~
   % cp ~/myaiida/aiida-vasp/examples/run_relax.py .
   % cp ~/myaiida/aiida-vasp/examples/auxiliary.py .

Usually a little modification of ``run_relax.py`` is necessary to run
this example, such as the queueing system job setting,

::

       options.resources = {'num_machines': 1,
                            'num_mpiprocs_per_machine': 16,
                            'tot_num_mpiprocs': 16}

Command options of ``run_relax.py`` are handled by the code written in
``auxiliary.py`` and the calculation is sent to AiiDA daemon by

::

   % python run_relax.py --potential-family PBE.54 vasp544mpi mycluster

The status of this calculation is checked by ``verdi calculation list``.

::

   % verdi calculation list
     PK  Creation    State           Type       Computer    Job state
   ----  ----------  --------------  ---------  ----------  --------------------
    684  15s ago     Waiting | None  vasp.vasp  mycluster   WITHSCHEDULER | None

   Total results: 1

   Info: last time an entry changed state: 2s ago (at 06:46:17 on 2018-11-15)

   % verdi calculation list
   PK    Creation    State    Type    Computer    Job state
   ----  ----------  -------  ------  ----------  -----------

   Total results: 0

   Info: last time an entry changed state: 0s ago (at 06:47:43 on 2018-11-15)

   % verdi calculation list -a -l 1
     PK  Creation    State         Type       Computer    Job state
   ----  ----------  ------------  ---------  ----------  ---------------
    684  2m ago      Finished | 0  vasp.vasp  mycluster   FINISHED | DONE

   Total results: 1

   Info: last time an entry changed state: 1m ago (at 06:47:43 on 2018-11-15)

Once the calculation is performed successively, the graph of created
nodes is watched quickly in a text form by

::

   % verdi node tree 684
   Property       Value
   -------------  ------------------------------------
   type           VaspCalculation
   pk             684
   uuid           c7f0005b-668d-4104-83ea-c43757a8ea8e
   label
   description
   ctime          2018-11-15 06:46:04.199212+00:00
   mtime          2018-11-15 06:47:43.035540+00:00
   process state  ProcessState.FINISHED
   exit status    0
   computer       [1] mycluster
   code           vasp544mpi

                           /-RemoteData [685]
                          |
   -- /VaspCalculation [684]-FolderData [686]
                          |
                           \-Dict [687]

Once a calculation could run successively, it is time to start trying
AiiDA tutorial (http://www.aiida.net/tutorials/) with AiiDA-VASP and
reading AiiDA documentation
(https://aiida-core.readthedocs.io/en/latest/). By this one
calculation, we can learn how to interact with our data using
``verdi`` command and python interactive shell (ipython invoked by
``verdi shell``). Although the amount of AiiDA documentation is large,
the most of them should be understood to draw workflow. Currently many
details of AiiDA are not yet written in the documentation. Therefore
we sometimes need to read the python docstrings to achieve our
workflows, for which the understanding of the fundamental concepts of
properties of AiiDA will be the key to find the entry points.
