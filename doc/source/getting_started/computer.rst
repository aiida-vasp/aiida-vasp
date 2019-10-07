.. _computer:

====================
8. Create a computer
====================

In order to execute any calculations, AiiDA needs a `computer`. This
can be a local computer or a cluster. Let us configure a
cluster and call it ``mycluster``. We will utilize SSH as the
transport (e.g. how `AiiDA`_ talks to the computer) and the `Torque`_
sheduler (`AiiDA`_ also supports the popular `Slurm`_ and `PBS`_).  In the
process you also need to specify the working directory on the cluster,
which is typically where you calculations are executed.
Typically, this is different from your home directory on your
cluster. Please consult your cluster's documentation to figure out the right
path.

Establish passwordless access to the computer
---------------------------------------------

Make sure you have passwordless access to ``mycluster``. We will not cover all details related
to this, but briefly explain the standard procedure. If you get stuck, consider Google as a
great resource to solve these issues. We can facility passwordless access by using keys.
You are in possession of a private key and a public key. The private key should never be
shared with anyone, not even another computer. However, the public key can (and should) be shared to
enable passwordless access.

#. We now assume you do not already have any SSH keys present. If that is the case you will be able
to figure this out yourself. Also the following strategy only works on standard Linux distributions.
First, create the key pair with ``ssh-keygen``:

::

   $ ssh-keygen
   Generating public/private rsa key pair.
   Enter file in which to save the key (/home/efl/.ssh/id_rsa):
   Enter passphrase (empty for no passphrase):
   Enter same passphrase again:
   Your identification has been saved in id_rsa.
   Your public key has been saved in id_rsa.pub.
   The key fingerprint is:
   SHA256:fU/YsA8uUp1gdCrPcGs12pWetNZpdLcT7DEic+MjTnA someuser@localhost
   The key's randomart image is:
   +---[RSA 2048]----+
   |                 |
   |        . .      |
   |       . o  . .  |
   |      o +o=E.O.=.|
   |       *SO=OO+*.=|
   |        *.O**= + |
   |       . ++oo o .|
   |        ....     |
   |         ..      |
   +----[SHA256]-----+

Please make sure you enter a reasonably long passphrase of your liking. You need to remember this to unlock the key pair in the future.

.. warning::

   Also consider that the standard RSA 2048 is not considered secure anymore. If ``mycluster``
   supports the ed25519 (elliptic curve) standard, please consider to use that instead by issuing
   ``ssh-keygen -t ed25519 -o -a 100``.

#. With this in place SSH into your computer to make sure this computer is a known host. Answer
   yes to the question. Exit and go back to your local computer.

#. We now need to put the public key (content of ``id_rsa.pub``) into the ``authorized_keys`` file
   that you find in the ``~/.ssh`` folder on ``mycluster``::

     ssh-copy-id -i ~/.ssh/id_rsa <username>@<mycluster.hostname>

   where ``<username>`` is your username and ``<mycluster.hostname>`` is the hostname of ``mycluster``.

#. Try to login to ``mycluster``. You should now (maybe after being asked to unlock the key pair)
   be able to access it without typing a password.

#. If you have problems, e.g. still have to enter the password please check the permissions
   of the ``~/.ssh/id_rsa`` file (local computer), ``~/.ssh/id_rsa.pub`` file (local computer),
   ``~/.ssh`` directory (both local and on server) and the ``~/.ssh/authorized_keys`` file (on server).
   Sometimes the directory does not exists at ``mycluster`` or has the wrong permissions. In fact,
   this is a rather common problem. You typically change the permissions of directories and files
   with the ``chmod`` command. Make sure the following permissions are set or present at ``mycluster`` ::

     chmod 700 ~/.ssh
     touch ~/.ssh/authorized_keys
     chmod 644 ~/.ssh/authorized_keys

   and::

     chmod 700 ~/.ssh
     chmod 600 ~/.ssh/id_rsa
     chmod 644 ~/.ssh/id_rsa.pub

   on your local computer and try again. If this does not work, talk to someone that is
   familiar with setting up this on your system environment.

We are now ready to add the computer.

Adding the actual computer to AiiDA
-----------------------------------

Let us now add the cluster computer to AiiDA by executing the following
commands::

   %/$ verdi computer setup
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

During the setup you will be asked for prepend and append text. Typically you leave these
empty. However, in case you would like to for instance configure several clusters and shift
between them you can enter specific directories here. Say for instance that on of the clusters need
a specific account number and the others no, then you typically enter this into the prepend text
of the given form required by your scheduler. When calling the calculation later, you would then
leave the account number empty such that the one you have defined in the prepend section is picked
up when loading the respective computers.

We are not entirely done, as we also need to configure the SSH
transport, which is done by::

   %/$ verdi computer configure ssh mycluster
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
administrator to get the correct details. Notice that we here demonstrated the use of
`?` to get more help and information at a given step. Finally, test that the computer ``mycluster``
works and is accessible from `AiiDA`_ by executing:

::

   %/$ verdi computer test mycluster
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

Great. Now the computer seems to work. Let us move on to the code.

.. _AiiDA: https://www.aiida.net
.. _Slurm: https://slurm.schedmd.com/
.. _Torque: https://www.adaptivecomputing.com/products/torque
.. _PBS: https://www.pbspro.org/
