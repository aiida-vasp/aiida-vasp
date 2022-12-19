Welcome to AiiDA-VASP's documentation!
======================================

AiiDA-VASP is a plug-in for the workflow management and data provenance tracking framework `AiiDA`_. It provides the classes `AiiDA`_ needs to run simulations using `VASP`_ (Vienna Ab initio Simulation Package). `VASP`_ is a program for atomic scale materials modelling, e.g. electronic structure calculations and quantum-mechanical molecular dynamics, from first principles. For detailed documentation on using `VASP`_ take a look in their `VASP wiki`_.

AiiDA-VASP is under active development, check out the `changelog`_.

.. warning::
   Please accept that the development of this plugin is a community effort and any error, bug or missing functionality that might appear is the responsibility of the individual user. If you detect something that needs to be improved we highly encourage to open an issue or even better, submit a pull requests, where you try to fix the issue yourself. You can do this on our Github repository.

.. note::
   We are currently looking for additional developers. If you are interested, please open an issue on our repository on Github.

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
   getting_started/install
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

   concepts/parameters
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

   developments/general
   developments/developments
   developments/running_tests
   changelog

.. toctree::
   :maxdepth: 1
   :caption: API reference
   :hidden:

   api/aiida_vasp
