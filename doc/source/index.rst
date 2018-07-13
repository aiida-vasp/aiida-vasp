.. AiiDA-VASP documentation master file, created by
   sphinx-quickstart on Mon Feb 12 09:24:45 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AiiDA-VASP's documentation!
======================================

AiiDA-VASP is a plug-in for the workflow management and data provenance tracking framework `AiiDA`_. It provides the classes `AiiDA`_ needs to run simulations using `VASP`_ (Vienna Ab initio Simulation Package). `VASP`_ is a program for atomic scale materials modelling, e.g. electronic structure calculations and quantum-mechanical molecular dynamics, from first principles.  For more info and a list of features look `here`_. For detailed documentation about using VASP take a look at the `documentation page`_ or the `wiki`_

AiiDA-VASP is under active development, check out the newest changes here: `changelog`_

.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at/
.. _here: https://www.vasp.at/index.php/about-vasp/59-about-vasp
.. _documentation page: https://www.vasp.at/index.php/documentation
.. _wiki: http://cms.mpi.univie.ac.at/wiki/index.php/The_VASP_Manual
.. _changelog: https://github.com/aiidateam/aiida-vasp/blob/develop/CHANGELOG.md


Getting started
===============

If you use python for other things chances are you have a working system to manage your virtual environments. Please note that AiiDA as well as this plug-in are tested with both ``conda`` as well as ``virtualenvwrapper``.

If you are already using AiiDA, simply activate the virtualenv you are using it in. Otherwise, set up a python virtualenv::

   $ pip install virtualenvwrapper
   $ mkvirtualenv --python=python2.7 aiida-vasp-env
   $ workon aiida-vasp-env

Or using conda::

   $ conda create -n aiida-vasp-env python=2
   $ source activate aiida-vasp-env

Install the plug-in using::

   $ pip install aiida-vasp
   $ reentry scan -r aiida  # (not necessary if using 'develop' branch of aiida)

This will automatically install the AiiDA python package(s) as well as any other dependencies of the plug-in and register all the plugin classes with AiiDA. Follow the steps in the `AiiDA documentation`_ to complete setting up AiiDA. Of course, if you had AiiDA already set up, you don't need to do that.

After setting up the database and profile and configuring the compute resources, you might want to run an example VASP calculation.

   $ (aiida-venv) git clone github.com/aiidateam/aiida-vasp
   $ (aiida-venv) python aiida-vasp/examples/run_vasp simple --import-from <POTCAR-path> <code> <computer>

Where ``<POTCAR-path>`` is the path to a set of POTCAR files (for example ``.../vasp_pot/potpaw_PBE``), ``<code>`` is the PK or name of the code you set up in AiiDA for running VASP, ``<computer>`` is the PK or name of the computer you set up in AiiDA for running VASP on.

.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/stable/

Running calculations
--------------------

 * Take a look at the file `example calc`_ for an example code on how to create and submit a VASP calculation from python code.
 * Take a look at the file `example workflow`_ for an example on how to do the same via an AiiDA WorkChain.
 * Read about running structure relaxations in the `How To section <howto/relax_wf/one-off>`

.. _example calc: https://github.com/aiidateam/aiida-vasp/blob/develop/examples/run_vasp.py
.. _example workflow: https://github.com/aiidateam/aiida-vasp/blob/develop/examples/run_base_wf.py

Managing potcar files
---------------------

AiiDA-VASP takes care of managing your POTCAR files, but because they are part of the VASP licence, you need to obtain them separately and make them available to AiiDA-VASP. You should have recieved a folder (``tar`` archive) containing multiple subfolders (``tar`` archives), each representing a set of POTCAR files intended to be used together. AiiDA-VASP allows you to upload only the sets (or even individual potentials) you require, and keep them grouped in so called "families".

The command line tools for these tasks are written as plugins to AiiDA, they can be called through AiiDA's ``verdi`` command like so::

   $ verdi data vasp-potcar --help
   Usage: verdi data vasp-potcar [OPTIONS] COMMAND [ARGS]...

      Top level command for handling VASP POTCAR files.

   Options:
     --help  Show this message and exit.
   
   Commands:
     exportfamily  Export a POTCAR family into a compressed tar...
     listfamilies  List available families of VASP potcar files.
     uploadfamily  Upload a family of VASP potcar files.

To make for example the PBE.54 family of POTCAR files available, use the ``uploadfamily`` command like so::

   $ verdi data vasp-potcar uploadfamily --path=vasp_pot/potpaw_PBE.54.tar --name=PBE.54 --description="PBE potentials for version 5.4"

Which will allow you to pass for example the following to the base workflow::

   $ inputs.potcar_family = Str('PBE.54')
   $ inputs.potcar_mapping = DataFactory('parameter')(dict={'In': 'In_d', 'As': 'As'})

Assuming you will run VASP on an InAs structure and wish to use the ``potpaw_PBE.54/In_d/POTCAR`` and the ``potpaw_Ppotpaw_PBE.54/As/POTCAR`` potentials.

More information about managing POTCAR files can be found here:

.. toctree::
   :maxdepth: 2

   howto/upload_potcars

Creating workflows
------------------
Read about how to write workflows by combining the building blocks provided as well as about the building blocks themselves.

.. toctree::
   :maxdepth: 3

   howto/write_workflows
   howto/use_relax_wf

More
====

.. The following may be partially outdated and is in the process of being brought up to date

.. toctree::
   :maxdepth: 4
   :caption: Contents:

   index_old



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
