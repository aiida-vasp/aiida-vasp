.. AiiDA-VASP documentation master file, created by
   sphinx-quickstart on Mon Feb 12 09:24:45 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AiiDA-VASP's documentation!
======================================

AiiDA-VASP is a plugin for the workflow management and data provenance tracking framework `AiiDA`_. It provides the classes `AiiDA`_ needs to run simulations using `VASP`_ (Vienna Ab initio Simulation Package). `VASP`_ is a program for atomic scale materials modelling, e.g. electronic structure calculations and quantum-mechanical molecular dynamics, from first principles.  For more info and a list of features look `here`_. For detailed documentaion about using VASP take a look at the `documentation page`_ or the `wiki`_


.. _AiiDA: https://www.aiida.net

.. _VASP: https://www.vasp.at/
.. _here: https://www.vasp.at/index.php/about-vasp/59-about-vasp
.. _documentation page: https://www.vasp.at/index.php/documentation
.. _wiki: http://cms.mpi.univie.ac.at/wiki/index.php/The_VASP_Manual

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

.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/stable/

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
