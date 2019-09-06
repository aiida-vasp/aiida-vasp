.. _getting_started:

=================
AiiDA VASP plugin
=================

.. image:: https://travis-ci.org/aiida-vasp/aiida-vasp.svg?branch=develop
   :target: https://travis-ci.org/aiida-vasp/aiida-vasp
			
.. image:: https://readthedocs.org/projects/aiida-vasp/badge/?version=latest
   :target: http://aiida-vasp.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
   
.. image:: https://coveralls.io/repos/github/aiida-vasp/aiida-vasp/badge.svg?branch=develop
   :target: https://coveralls.io/github/aiida-vasp/aiida-vasp?branch=develop
      
This is a plugin to `AiiDA`_ to run calculations with the ab-initio program `VASP`_.

Please have a look at the `Documentation <https://aiida-vasp.readthedocs.org/en/latest>`_ for additional details.

Starting to use the plugin
--------------------------

1. If you are already using `AiiDA`_, simply activate the virtual environment associated with it, here assumed to be located in ``~/env/aiida-vasp``::
     
   $ source ~/env/aiida-vasp/bin/activate

And skip to step 4 and continue from there.

2. Otherwise, set up a new virtual environment::

   $ python -m venv ~/env/aiida-vasp

3. Enable the newly installed virtual environment::

   $ source ~/env/aiida-vasp/bin/activate

4. Install the plugin::

   $ (aiida-vasp) pip install aiida-vasp

5. Update the entry points that `AiiDA`_ are using::

   $ (aiida-vasp) reentry scan -r aiida

This will automatically install the `AiiDA`_ python package(s) as well as any other dependencies of the plugin and register all the plugin classes with `AiiDA`_. Follow the steps in the `AiiDA documentation`_ to complete setting up `AiiDA`_ in case you have not yet completed these steps. This involves setting up a `profile`, a `computer` and a `code`.

After setting up these, you might want to run an example `VASP_` calculation to check that everything works as expected. In order to do so, we need to fetch the examples from the GitHub repository. We believe it is good practice to put things related to a specific virtual environment in its associated folders.

6. Download the full AiiDA-VASP repository, which includes the examples::

   $ (aiida-vasp) cd ~/env/aiida-vasp
   $ (aiida-vasp) git clone https://github.com/aiida-vasp/aiida-vasp.git

The examples are located inside the ``examples`` directory. The most simple example ``run_vasp_lean`` is what you need for testing at this point.

But before we can execute a calculation, we need to upload the potentials, or the POTCAR files. We know you are keen on getting the first calculation done, but how the plugin handles the potentials are important to understand as they are covered by the license. Please consult :ref:`potentials`.

7. We now assume that you are in the root folder of a offically supported potential family bundle supplied by the `VASP`_ team called ``paw_gga_pbe_54`` located in the folder ``~/potentials`` of the user. This potential family can now be uploaded using::
  
   $ cd ~/potentials/paw_gga_pbe_54
   $ (aiida-vasp) verdi data vasp-potcar uploadfamily --name=PBE.54 --description="The PBE 5.4 potentials"


The plugin will then traverse the folder you are in and look for POTCAR files. If you want you could also pass an argument that points to the folder, instead of changing the directory, e.g::
  
   $ (aiida-vasp) verdi data vasp-potcar uploadfamily --path=~/potentials/paw_gga_pbe_54 --name=PBE.54 --description="The PBE 5.4 potentials

Or you could also just use the ``tar`` file directly::

   $ (aiida-vasp) verdi data vasp-potcar uploadfamily --path=~/paw_gga_pbe_54.tar --name=PBE.54 --description="The PBE 5.4 potentials
   
This will probably take some time. Give it time to finish. You will get a warning that the formatting is off for one of the POTCAR_H files, but you can ignore that. Additional details regarding the handling of the POTCAR files can be found here :ref:`potentials`.

8. Finally we are ready to launch a test calculation. Change to the ``examples`` directory in the plugin's root directory and run the command::
     
   $ (aiida-vasp) cd ~/env/aiida-vasp/aiida-vasp/examples
   $ (aiida-vasp) python run_vasp_lean.py <code> <computer> --potential-family=PBE.54

Where ``<code>`` is the PK or name of the `code` you set up in `AiiDA`_ for running `VASP`_ and ``<computer>`` is the PK or name of the `computer` you set up in AiiDA for running `VASP`_ on. The attribute ``potential-family`` is the name you used when uploading the potential family.

In the ``examples`` folder you will find additional example files that executes certain sets of workchains, which again executes `VASP`_ accordingly when needed. Currently, the following is bundled with the distribution:

 * ``run_vasp_lean`` - a clean execution of `VASP`_ with minimal functionality
 * ``run_relax`` - runs structure relaxations on top of ``run_relax``
 * ``run_converge`` - runs k-point and plane wave cutoff convergence tests on top of ``run_relax``
 * ``run_bands`` - an example of how to use the ``master`` workchain to perfom a more composed set of calculation with dependencies between each step (in this case ``CHGCAR`` among other things)

The idea is that one should always run ``run_converge`` for all calculations. That serves as an entry point to a `VASP`_ calculation. As you can see in the documentation for the convergence and relaxation workchain, the actual convergence calculation can be skipped if the necessary parameters are given.

Please consult the navigation bar for additional details on how to use the plugin.

.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/latest/
