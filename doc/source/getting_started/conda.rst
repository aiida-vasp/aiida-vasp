.. _conda:

============================
Setup in a Conda environment
============================

AiiDA-core and AiiDA-VASP may be installed on your home
directory. Because AiiDA core and AiiDA-VASP are still under the
acitive development, they are updated frequently and you may be needed
to test across different versions. It is a good idea to manage your
AiiDA systems of different versions and/or settings separately. This
becomes possible to use a virtual environment in Python or the conda
environment functionality.

In this guide, the conda environment is used to explain to establish
creating your AiiDA system. The reason to use Conda is that it is easy
to manage and install PostgreSQL as a non-privileged user being a
requirement of AiiDA core (hereafter AiiDA). Here RabbitMQ, that is
also a requirement of AiiDA, is assumed already installed on the
system.

In the following, the preparation of AiiDA and AiiDA-VASP on a Linux
system is explained in a step-by-step fashion. We assume Conda is
installed and operational.

Preparing the Conda environment
-------------------------------

Each installation of AiiDA (and its settings) is isolated by first
creating a containing directory (in the following we call this
``myaiida`` placed in the root home folder of the active user). In it,
AiiDA, AiiDA-VASP plugin, and PostgreSQL database are installed:

::

   % cd ~
   % mkdir myaiida
   % cd myaiida
   % conda create -n myaiida python=3
   % source activate myaiida
   % conda install gcc_linux-64 gxx_linux-64 ipython

``conda activate myaiida`` enables the Conda environment ``myaiida``.
All the settings, installs or running are always done in this conda
environment. It is mandatory to activate this every time logging to
this computer. Please also remember to activate this everytime you
need to make changes to the respective environments. You can of course
created multiple Conda environments depending on your use case or versions
you would like to have installed.


Setting up the PostgreSQL database
----------------------------------
Next we install and setup PostgreSQL. Below, ``mypassword`` is used as
the PostgreSQL user password for a database user ``aiida``. For the
time being, we assume that the database user is not the same or
correlated with the system users. Please also consider that
``mypassword`` is not hidden and stored in clear text. It is thus
recommended to use a non-critical password that you do not use for
something else. Let us install and setup the first database called
``aiidadb`` by issuing the following commands:

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
invoked by the `psql` command. If a user wants, it is fully possible
to have sevaral AiiDA databases running (also using the same user and
password). In AiiDA it is possible to change which database to
use. Keep that in mind, in case a use case should present itself to
you, where this makes sense.

Install AiiDA and setup a profile
---------------------------------

We are now ready to install AiiDA. Here we use AiiDA release of
`v1.0.0b6`. Until AiiDA makes an offical `v1.x` release we recommend
to clone the Git repository and use the developer install option for
``pip``. When doing this, when you later issue ``git pull`` in the
AiiDA directory, your installation is automatically updated (unless
there has been any big changes, which are not likely to happen at this
point).  We can install AiiDA by issuing the following commands (make
sure you are still in the root of the ``myaiida`` directory, which you
should be, unless you have diverted from this guide):

::

   % git clone https://github.com/aiidateam/aiida-core.git
   % cd aiida_core
   % git checkout v1.0.0b6
   % pip install -e .

Now we create the AiiDA setup with the profile name ``myaiida``
(i.e. the same name as used for the Conda environment folder) as
follows (replace whatever suits your needs)::

   % export AIIDA_PATH=~/myaiida
   % verdi setup
   Info: enter "?" for help
   Profile name: myaiida
   User email: mymail@address.com
   First name: MyFirstName
   Last name: MyLastName
   Institution: My university
   Password [***]:
   Database engine (postgresql_psycopg2) [postgresql_psycopg2]:
   Database backend (django, sqlalchemy) [django]:
   Database hostname [localhost]:
   Database port [5432]:
   Database name: aiidadb
   Database username: aiida
   Database password: mypassword
   Repository directory [/home/username/myaiida/.aiida/repository/migration]:

The command ``verdi`` controlls AiiDA and gives access points from the
command line.  Of course you can also use AiiDA directly from Python,
but sometimes ``verdi`` is easier and faster to use. To use AiiDA
``verdi`` command with tab completion, the file
``$CONDA_PREFIX/etc/conda/activate.d/env_vars.sh`` needs to be created
with the following content::

   #!/bin/sh
   export AIIDA_PATH=~/myaiida
   eval "$(_VERDI_COMPLETE=source-bash verdi)"
   export HOST=`hostname`

and ``$CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh`` containing::

   #!/bin/sh
   export HOST=`hostname`

In order to enable the recently placed files, we need to reactivate
the Conda environment by issuing the following:

::

   % conda deactivate
   % conda activate myaiida

AiiDA relies on a daemon that monitors and controlls your
calculations. You should now be able to start the AiiDA daemon by
issuing:

::

   % verdi daemon start

Sometimes this command fails. If it does, run ``verdi setup myaiida`` and
``verdi daemon start`` again.


Install AiiDA-VASP plugin
-------------------------

The AiiDA-VASP plugin does not need to be configured, or a profile created. It is
simply just an install, which is done as follows::

   % pip install aiida-vasp
   % reentry scan -r aiida

Setup up a computer in AiiDA
-----------------------------

In order to execute any calculations, AiiDA needs a ``computer``. This
can be a local computer, cluster, super computer. Let us configure a
cluster and call it ``mycluster``. We will utilize SSH as the
transport (e.g. how AiiDA talks to the computer) and the Torque
sheduler (AiiDA also supports the popular Slurm and PBS).  In the
process you also need to specify the working directory on the cluster,
which is typically where you calculations are executed on the
cluster. Typically, this is different from your home directory on your
cluster. Remember you can enter `?` to get help at any point. Let us
now add the cluster computer to AiiDA by executing the following
commands:

::

   % verdi computer setup
   Info: enter "?" for help
   Computer label: mycluster
   Hostname: mycluster
   Description []:
   Enable the computer? [True]:
   Transport plugin: ssh
   Scheduler plugin: torque
   Shebang line (first line of each script, starting with #!) [#!/bin/bash]:
   Work directory on the computer [/scratch/{username}/aiida/]: /home/username/aiida/
   Mpirun command [mpirun -np {tot_num_mpiprocs}]:
   Success: Computer<1> mycluster created
   Info: Note: before the computer can be used, it has to be configured with the command:
   Info:   verdi computer configure ssh mycluster

We are not entirely done, as we also need to configure the SSH
transport, which is done by:

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

These parameters should be self explanatory. In case of doubt please contant your IT
administrator to get the correct details. Make sure that the active system user have
keyless access to the cluster. Finally, test that the computer ``mycluster``
works and is accessible from AiiDA by

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


Setup a code in AiiDA
---------------------

Now we need to add the code (in this case VASP) to AiiDA.  The
subcommand ``code`` describes by which code our calculations run. The
code has to be installed on the location of ``computer``, i.e., if it
is a computer cluster, the code has to be installed properly to run
there. The setup is done as follows::

   % verdi code setup
   Info: enter "?" for help
   Label: vasp544mpi
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
   Remote absolute path: /usr/local/calc/vasp/vasp544mpi
   Success: Code<1> vasp544mpi@mycluster created

We can check if the code is now listed by issuing:

::

   % verdi code list
   # List of configured codes:
   # (use 'verdi code show CODEID' to see the details)
   * pk 1 - vasp544mpi@mycluster
   # No codes found matching the specified criteria.

And look at its details. These commands are also available for the computers.

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


Run an AiiDA-VASP calculation
-----------------------------

AiiDA relies mainly on the concept of ``workchain`` which is a
composition of a setup and teardown of a calculation (or calls to
other ``workchain``).  A ``workchain`` can be composed into one or
multiple `workflows`. A small amount of basic ``workchain``'s are
included in AiiDA-VASP. Users are encouraged to develop new, or
complementig ``workchain``'s and submitting them to the repository to
increase the efficiency of all VASP users.

An example of a ``workchain`` calculation, which performs relaxation,
is copied from the ``example`` directory.

::

   % cd ~
   % mkdir run_example && cd run_example
   % cp ~/myaiida/aiida-vasp/examples/run_relax.py .
   % cp ~/myaiida/aiida-vasp/examples/auxiliary.py .

Usually a little modification of ``run_relax.py`` is necessary to run
this example, such as the queueing system job setting:

::

       options.resources = {'num_machines': 1,
                            'num_mpiprocs_per_machine': 16,
                            'tot_num_mpiprocs': 16}

maybe also setting the ``qos`` or the ``account`` etc., see the
available parameters at `AiiDA documentation
<https://aiida.readthedocs.io/projects/aiida-core/en/latest/scheduler/index.html>`_.

Command options of ``run_relax.py`` are handled by the code written in
``auxiliary.py`` and the calculation is sent to AiiDA daemon by executing:

::

   % python run_relax.py --potential-family PBE.54 vasp544mpi mycluster

We thus tell AiiDA, which potential family, code and computer to use.
The status of the execution of this ``workchain`` can be checked by
``verdi process list``.

::

   % verdi report list
     PK  Created    Process label      Process State    Process status
   ----  ---------  -----------------  ---------------  ----------------------------------
   6637  15s ago    RelaxWorkChain      Waiting         Waiting for child processes: 6640
   6640  10s ago    VerifyWorkChain     Waiting         Waiting for child processes: 6641
   6641  7s ago     VaspWorkChain       Waiting         Waiting for child processes: 6642
   6642  2s ago     VaspCalculation     Waiting         Waiting for transport task: upload

When executing ``run_relax.py``, in fact, three ``workchain`` are
executed. This is typically how you build workflows. In this case,
only ``VaspWorkChain`` calls a ``VaspCalculation`` process, which
again is responsible for calling VASP itself. When the execution is
complete, the graph can be created and inspected.

::



Once the example calculation above executed successively, it is time
to start trying AiiDA tutorial (http://www.aiida.net/tutorials/) with
AiiDA-VASP and reading AiiDA documentation
(https://aiida-core.readthedocs.io/en/latest/). By using this example
calculation, we can learn how to interact with our data using
``verdi`` command and python interactive shell (ipython invoked by
``verdi shell``). Although the amount of AiiDA documentation is large,
it should be understood from a viewpoint of to designing workflows and
managing data. That is after all the main purpose of AiiDA. Currently
many details of AiiDA are not yet documentated.
