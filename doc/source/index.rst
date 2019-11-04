 .. AiiDA-VASP documentation master file, created by
   sphinx-quickstart on Mon Feb 12 09:24:45 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AiiDA-VASP's documentation!
======================================

AiiDA-VASP is a plug-in for the workflow management and data provenance tracking framework `AiiDA`_. It provides the classes `AiiDA`_ needs to run simulations using `VASP`_ (Vienna Ab initio Simulation Package). `VASP`_ is a program for atomic scale materials modelling, e.g. electronic structure calculations and quantum-mechanical molecular dynamics, from first principles. For detailed documentation on using `VASP`_ take a look in their `VASP wiki`_.

AiiDA-VASP is under active development, check out the `changelog`_.

The plugin is only tested with VASP 5.4.4 and is compatible with Python>=3.5.

.. warning:: If you export any calculation done with the plugin the POTCAR file will be in the export file. Please remember that
   this file is protected by license and make sure all your actions in this process does not void this license.
   We take no responsibility and making sure the license if obeyed at all times is a user responsibility. We are currently
   working on a solution where it is going to be possible to export data that does not contain the POTCAR file (but
   it will re-link itself to the right POTCAR content on import, assuming the exact same file is present in the database).

General comments about the user documentation and support
=========================================================

.. _convention:

We will try to keep the documentation updated for both the standard Python and the `Conda`_ virtual environment. In order to make this as clear as possible, we have opted for the following syntax::

  $ commands in a standard Linux console, used together with the standard Python virtual environment

  % commands in the Conda console, used together with the Conda virtual environment

  %/$ same commands in the Conda console or standard Linux console

  # commands executed in the Postgresql interactive terminal

As such, whenever needed, we will try to keep parallel documentation on this. If you are on a Conda
platform, please follow the commands prepended with ``%``. Similarly you follow the commands with
a ``$`` prepended if you are using the standard Python virtual enviroments. In case the commands are
the same for both cases, ``%/$`` is prepended to the commands. The commands executed in
the Postgresql interactive terminal (``psql``) should be similar between the two.

.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _VASP wiki: https://cms.mpi.univie.ac.at/wiki/index.php
.. _changelog: https://github.com/aiida-vasp/aiida-vasp/blob/develop/CHANGELOG.rst
.. _Conda: https://docs.conda.io/en/latest/

.. toctree::
   :maxdepth: 1
   :caption: Getting started
   :hidden:

   getting_started/general
   getting_started/python
   getting_started/conda
   getting_started/rabbitmq
   getting_started/database
   getting_started/install
   getting_started/interact
   getting_started/profile
   getting_started/computer
   getting_started/code
   getting_started/potentials
   getting_started/test_run

.. toctree::
   :maxdepth: 1
   :caption: Tutorials
   :hidden:

   tutorials/fcc_si_step1
   tutorials/fcc_si_step2
   tutorials/fcc_si_step3
   tutorials/fcc_si_dos
   tutorials/run_vasp_builder
   tutorials/bulk_modulus
   tutorials/bulk_modulus_2
   tutorials/band_structure_electrons
   tutorials/a_daily_life

.. toctree::
   :maxdepth: 1
   :caption: Concepts
   :hidden:

   concepts/potentials
   concepts/calculations
   concepts/workchains
   concepts/workflows
   concepts/parsing
   concepts/immigrations

.. toctree::
   :maxdepth: 1
   :caption: Workflows
   :hidden:

   workflows/using_workflows
   workflows/designing_workflows

.. toctree::
   :maxdepth: 1
   :caption: Calculations
   :hidden:

   calculations/vasp
   calculations/wannier
   calculations/immigrator

.. toctree::
   :maxdepth: 1
   :caption: Workchains
   :hidden:

   workchains/vasp
   workchains/verify
   workchains/relax
   workchains/converge
   workchains/bands
   workchains/master
   workchains/writing_workchains

.. toctree::
   :maxdepth: 1
   :caption: Developments
   :hidden:

   developments/developments
   developments/running_tests
   changelog

.. toctree::
   :maxdepth: 1
   :caption: API reference
   :hidden:

   api/aiida_vasp
