.. _getting_started:

=================
AiiDA VASP plugin
=================

.. |version| image:: https://img.shields.io/pypi/v/aiida-vasp
   :target: https://pypi.org/project/aiida-vasp/
   :alt: Stable version

.. |status| image:: https://img.shields.io/pypi/status/aiida-vasp
   :target: https://pypi.org/project/aiida-vasp/
   :alt: PyPI - Status

.. |versions| image:: https://img.shields.io/pypi/pyversions/aiida-vasp
   :target: https://pypi.org/project/aiida-vasp/
   :alt: Supported Python versions

.. |build| image:: https://github.com/aiida-vasp/aiida-vasp/workflows/aiida-vasp/badge.svg
   :target: https://github.com/aiida-vasp/aiida-vasp/action
   :alt: Build status

.. |coverage| image:: https://codecov.io/gh/espenfl/aiida-vasp/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/espenfl/aiida-vasp
   :alt: Test coverage

.. |doc| image:: https://readthedocs.org/projects/aiida-vasp/badge/?version=latest
   :target: http://aiida-vasp.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |downloads| image:: https://img.shields.io/pypi/dm/aiida-vasp
   :target: https://pypi.org/project/aiida-vasp/
   :alt: PyPI - Downloads/month

.. |commits| image:: https://img.shields.io/github/commit-activity/m/aiida-vasp/aiida-vasp
   :target: https://github.com/aiida-vasp/aiida-vasp/commits/develop
   :alt: GitHub commit activity

+---------+-------------------------------+
| Release | |version| |status| |versions| |
+---------+-------------------------------+
| Build   | |build| |coverage| |doc|      |
+---------+-------------------------------+
| Stats   | |downloads| |commits|         |
+---------+-------------------------------+


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

   $ (aiida-vasp) pip install aiida-vasp

#. Update the entry points that `AiiDA`_ are using::

   $ (aiida-vasp) reentry scan -r aiida

This will automatically install the `AiiDA`_ python package(s) as well as any other dependencies of the plugin and register all the plugin classes with `AiiDA`_.

Please consider that `AiiDA`_ have prerequisite that needs to be installed and ensured working. The steps above will not take care of this for you. Please consult `AiiDA prerequisites`_ and follow the instructions therein.


Upgrading from AiiDA 1.x
------------------------

At the moment, only the development version of the plugin supports ``aiida-core >= 2.0.1``. 
If you are upgrading an existing AiiDA 1.x installation. Please upgrade ``aiida-core`` first and reinstall the plugin using::


    $ pip install git+https://github.com/aiida-vasp/aiida-vasp.git@develop#egg=aiida-vasp


This ensures the installation of the latest development version and registers the entrypoints. 
You can verify the latter with::

    $ verdi plugin list aiida.groups

and there should be entries of ``vasp.potcar``. The ``aiida.workflows`` and ``aiida.calculations`` entrypoints can be checked in a similiar way. 



.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/latest/
.. _AiiDA-VASP documentation: https://aiida-vasp.readthedocs.io/en/latest/
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _AiiDA prerequisites: https://aiida-core.readthedocs.io/en/latest/install/prerequisites.html
