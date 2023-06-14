Welcome to AiiDA-VASP's documentation!
======================================

AiiDA-VASP is a plug-in for the workflow management and data provenance tracking framework `AiiDA`_. It provides the classes `AiiDA`_ needs to run simulations using `VASP`_ (Vienna Ab initio Simulation Package). `VASP`_ is a program for atomic scale materials modelling, e.g. electronic structure calculations and quantum-mechanical molecular dynamics, from first principles. For detailed documentation on using `VASP`_ take a look in their `VASP wiki`_.

AiiDA-VASP is under active development, check out the `changelog`_.

.. warning::
   Please accept that the development of this plugin is a community effort and any error, bug or missing functionality that might appear is the responsibility of the individual user. If you detect something that needs to be improved we highly encourage to `open an issue`_ or even better, `submit a pull request`_, where you try to fix the issue yourself. You can do this on our Github repository.

.. note::
   We are currently looking for additional developers. If you are interested, please open an issue on our repository on Github.

.. note::
   Please consider to `open an issue`_ instead of sending an email to the AiiDA mailing list if the issue is related to this plugin.

Also, please consider that `AiiDA-VASP`_ is no substitute for not knowing how to use VASP. If in doubt at any point of the `VASP`_ parts of the
training material for `AiiDA-VASP`_ or when using `AiiDA-VASP`_, consult for instance `VASP lectures`_, `VASP tutorials`_, `VASP howtos`_,
`VASP tutorials using notebooks`_ or `VASP videos`_ or ask experienced `VASP`_ users.

.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _VASP wiki: https://cms.mpi.univie.ac.at/wiki/index.php
.. _changelog: https://github.com/aiida-vasp/aiida-vasp/blob/develop/CHANGELOG.rst
.. _Conda: https://docs.conda.io/en/latest/
.. _VASP lectures: https://www.vasp.at/wiki/index.php/Lectures_and_presentations
.. _VASP tutorials: https://www.vasp.at/wiki/index.php/Category:Tutorials
.. _VASP howtos: https://www.vasp.at/wiki/index.php/Category:Howto
.. _VASP tutorials using notebooks: https://www.vasp.at/tutorials/latest/
.. _VASP videos: https://www.youtube.com/channel/UCBATkNZ7pkAXU9tx7GVhlaw
.. _open an issue: https://github.com/aiida-vasp/aiida-vasp/issues
.. _submit a pull request: https://github.com/aiida-vasp/aiida-vasp/pull
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp

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
   :caption: Basic tutorials
   :hidden:

   tutorials/fcc_si_step1
   tutorials/fcc_si_step2
   tutorials/fcc_si_step3
   tutorials/fcc_si_step4
   tutorials/fcc_si_dos
   tutorials/interacting_data.rst

.. toctree::
   :maxdepth: 1
   :caption: Specific tutorials
   :hidden:

   tutorials/run_vasp_builder
   tutorials/bulk_modulus

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

.. toctree::
   :maxdepth: 1
   :caption: Calculations
   :hidden:

   calculations/vasp
   calculations/wannier

.. toctree::
   :maxdepth: 1
   :caption: Workchains
   :hidden:

   workchains/vasp
   workchains/relax
   workchains/converge
   workchains/bands
   workchains/writing_workchains

.. toctree::
   :maxdepth: 1
   :caption: Developments
   :hidden:

   developments/general
   developments/contributions
   developments/documentation
   developments/tests
   developments/vasp
   changelog

.. toctree::
   :maxdepth: 1
   :caption: API reference
   :hidden:

   apidoc/aiida_vasp
