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

Please have a look at the `AiiDA-VASP documentation`_ for instructions on how to install and use the plugin.

Installing the plugin
---------------------

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

This will automatically install the `AiiDA`_ python package(s) as well as any other dependencies of the plugin and register all the plugin classes with `AiiDA`_. Please consider that `AiiDA`_ have prerequisite that needs to be installed and ensured working. The steps above will not take care of this for you. Please consult `AiiDA prerequisites`_ and follow the instructions therein. Typically, on a standard Linux installation, like a modern Ubuntu, Fedora or similar distribution, the only thing you need to install is RabbitMQ and install and setup a Postgresql database. Please consider to work on a bleeding distribution of Fedora or Ubuntu and not CentOS or LTS versions. This also works fine, but requires more knowledge of the users to update system packages. After this is done, a `profile`, a `computer` and a `code` needs to be configured. The rest of the documentation gives information how to perform these steps (and only on Fedora 30 if there are system specific commands issued). If at any point you need more information, please consult the general `AiiDA-VASP documentation`_. In case you have already configured a `profile` and a `computer` please skip to the step which creates a `code`. It is highly unlikely that you need help with the installation of AiiDA-VASP and already have a VASP code configured in AiiDA. In case that is indeed so, please skip this step and launch your first test calculation and then continue on with the tutorial steps to learn more about the plugin.

.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/latest/
.. _AiiDA-VASP documentation: https://aiida-vasp.readthedocs.io/en/latest/
.. _AiiDA prerequisites: https://aiida-core.readthedocs.io/en/latest/install/prerequisites.html
