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

..
  Keep this comment as it is used for including these steps in the install section of the docs.
  It includes everything past the next line.
  Start installation description

#. If you are already using `AiiDA`_, simply activate the virtual environment associated with it, here assumed to be located in ``~/env/aiida-vasp``::

   $ source ~/env/aiida-vasp/bin/activate

#. Otherwise, set up a new virtual environment::

   $ python -m venv ~/env/aiida-vasp

#. And then enable the newly installed virtual environment::

   $ source ~/env/aiida-vasp/bin/activate

#. Install the `AiiDA-VASP`_ plugin (and `AiiDA`_ if that is not already installed)::

   $ (aiida-vasp) pip install aiida-vasp

  If you need to install the compatibility release of `AiiDA-VASP`_ which works with `AiiDA`_ 1.6.4 you should instead install the plugin
  using ``pip install aiida-vasp=2.2``, but this is not recommended and only mentioned for legacy support. For the legacy version you
  also most likely have to run ``reentry scan -r aiida`` after installing the plugin.

This will automatically install the `AiiDA`_ python package(s) as well as any other dependencies of the plugin and register all the plugin classes with `AiiDA`_.

Please consider that `AiiDA`_ have prerequisite that needs to be installed and ensured working. The steps above will not take care of this for you. Please consult `AiiDA prerequisites`_ and follow the instructions therein.

..
  End installation description

Support
-------

..
  Start support description

The development, maintenance and use of this plugin is considered a community effort. In order to facilitate for the community to contribute,
we have established a `space on Matrix`_ that users can use to communicate. We encourage users to help each other. In addition,
the development team is present in the space and users are free to ask.
First consult the documentation of both `AiiDA-VASP documentation`_ and `AiiDA documentation`_ and also consider that the developers are
not paid for this work. Please respect potential lead times in getting answers and be polite.

..
  End support description

.. _space on Matrix: https://matrix.to/#/#aiida-vasp:matrix.org
.. _AiiDA-VASP documentation: https://aiida-vasp.readthedocs.io/en/latest/
.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/latest/
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _AiiDA prerequisites: https://aiida-core.readthedocs.io/en/latest/install/prerequisites.html
