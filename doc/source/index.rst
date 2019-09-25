.. AiiDA-VASP documentation master file, created by
   sphinx-quickstart on Mon Feb 12 09:24:45 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AiiDA-VASP's documentation!
======================================

AiiDA-VASP is a plug-in for the workflow management and data provenance tracking framework `AiiDA`_. It provides the classes `AiiDA`_ needs to run simulations using `VASP`_ (Vienna Ab initio Simulation Package). `VASP`_ is a program for atomic scale materials modelling, e.g. electronic structure calculations and quantum-mechanical molecular dynamics, from first principles. For detailed documentation on using `VASP`_ take a look in their `VASP_wiki`_.

AiiDA-VASP is under active development, check out the :changelog:

General comments about the user documentation and support
=========================================================

If you already use Python chances are you have a working system to manage your virtual environments. Please note that `AiiDA`_ as well as this plugin are tested with both ``virtualenv``, ``conda`` and ``virtualenvwrapper``. Since Python 2 is soon to be deprecated the documentation are written assuming Python 3.6 or later is present and working. However, `AiiDA`_ and AiiDA-VASP are compatible with Python 2 in case you prefer. If that is the case you should be able to translate the procedures that are different to Python 2. In the documentation we will assume you are working with ``virtualenv``. which is also bundled with Python 3.

Tested with VASP 5.4.4.

.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _VASP_wiki: https://cms.mpi.univie.ac.at/wiki/index.php
.. _changelog: https://github.com/aiida-vasp/aiida-vasp/blob/develop/CHANGELOG.md


.. toctree::
   :maxdepth: 1
   :caption: Getting started
   :hidden:

   getting_started/index
   getting_started/conda

.. toctree::
   :maxdepth: 1
   :caption: Tutorials
   :hidden:

   tutorials/oslo
   tutorials/a_daily_life
   tutorials/run_vasp
   tutorials/relaxations
   tutorials/convergence_tests
   tutorials/band_structure_electrons
   tutorials/band_structure_phonons
   tutorials/bulk_modulus
   tutorials/bulk_modulus_2

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
   changelog

.. toctree::
   :maxdepth: 1
   :caption: API reference
   :hidden:

   api/aiida_vasp
