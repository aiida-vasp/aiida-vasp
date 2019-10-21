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

#. If you are already using `AiiDA`_, simply activate the virtual environment associated with it, here assumed to be located in ``~/env/aiida-vasp``::

   $ source ~/env/aiida-vasp/bin/activate

#. Otherwise, set up a new virtual environment::

   $ python -m venv ~/env/aiida-vasp

#. And then enable the newly installed virtual environment::

   $ source ~/env/aiida-vasp/bin/activate

#. Install the `AiiDA-VASP`_ plugin (and `AiiDA`_ if that is not already installed)::

   $ (aiida-vasp) pip install aiida-vasp --pre

#. Update the entry points that `AiiDA`_ are using::

   $ (aiida-vasp) reentry scan -r aiida

This will automatically install the `AiiDA`_ python package(s) as well as any other dependencies of the plugin and register all the plugin classes with `AiiDA`_.

Please consider that `AiiDA`_ have prerequisite that needs to be installed and ensured working. The steps above will not take care of this for you. Please consult `AiiDA prerequisites`_ and follow the instructions therein.

.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/latest/
.. _AiiDA-VASP documentation: https://aiida-vasp.readthedocs.io/en/latest/
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _AiiDA prerequisites: https://aiida-core.readthedocs.io/en/latest/install/prerequisites.html
