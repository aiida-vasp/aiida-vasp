.. AiiDA-VASP documentation master file, created by
   sphinx-quickstart on Mon Feb 12 09:24:45 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AiiDA-VASP's documentation!
======================================

AiiDA-VASP is a plugin for the workflow management and data provenance tracking framework `AiiDA`_. It provides the classes `AiiDA`_ needs to run simulations using `VASP`_ (Vienna Ab initio Simulation Package). `VASP`_ is a program for atomic scale materials modelling, e.g. electronic structure calculations and quantum-mechanical molecular dynamics, from first principles.  For more info and a list of features look `here`_. For detailed documentaion about using VASP take a look at the `documentation page`_ or the `wiki`_

AiiDA-VASP is under active development, check out the newest changes here: `changelog_`

.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at/
.. _here: https://www.vasp.at/index.php/about-vasp/59-about-vasp
.. _documentation page: https://www.vasp.at/index.php/documentation
.. _wiki: http://cms.mpi.univie.ac.at/wiki/index.php/The_VASP_Manual
.. _changelog: https://github.com/DropD/aiida-vasp/blob/develop/CHANGELOG.md


Getting started
===============

If you use python for other things chances are you have a working system to manage your virtual environments. Please note that AiiDA as well as this plugin are tested with both conda as well as virtualenvwrapper.

If you are already using AiiDA, simply activate the virtualenv you are using it in. Otherwise, set up a python virtualenv::

   $ pip install virtualenvwrapper
   $ mkvirtualenv --python=python2.7 aiida
   $ workon aiida-vasp

Install the plugin using::

   $ pip install git+https://github.com/DropD/aiida-vasp/tree/master#egg=aiida-vasp-0.1.0

This will automatically install the AiiDA python package(s) as well as any other dependencies of the plugin. Follow the steps in the `AiiDA documentation`_ to complete setting up AiiDA. Of course, if you had aiida already set up, you don't need to do that.

After setting up the database and profile and configuring the compute resources, you might want to run an example VASP calculation.

   $ (aiida-venv) git clone github.com/DropD/aiida-vasp
   $ (aiida-venv) python aiida-vasp/examples/run_vasp simple --import-from <POTCAR-path> <code> <computer>

Where ``<POTCAR-path>`` is the path to a set of POTCAR files (for example ``.../vasp_pot/potpaw_PBE``), ``<code>`` is the PK or name of the code you set up in AiiDA for running VASP, ``<computer>`` is the pk or name of the computer you set up in AiiDA for running VASP on.

Take a look at the file ``aiida-vasp/examples/run_vasp.py`` for example code on how to programmatically create and submit a VASP calculation.

.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/stable/


More
====

The following is documentation for a slightly out of date version, which was written to work with AiiDA up to version 0.7.

.. toctree::
   :maxdepth: 4
   :caption: Contents:

   index_old



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
