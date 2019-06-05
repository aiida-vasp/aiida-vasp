=================
AiiDA Vasp Plugin
=================

.. image:: https://travis-ci.org/aiida-vasp/aiida-vasp.svg?branch=develop
    :target: https://travis-ci.org/aiida-vasp/aiida-vasp

.. image:: https://readthedocs.org/projects/aiida-vasp/badge/?version=latest
   :target: http://aiida-vasp.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://coveralls.io/repos/github/aiida-vasp/aiida-vasp/badge.svg?branch=develop
   :target: https://coveralls.io/github/aiida-vasp/aiida-vasp?branch=develop

This is a plugin to `AiiDA <www.aiida.net/?page_id=264>`_ to run calculations with the ab-initio program `VASP <https://www.vasp.at/>`_.

* `Documentation <https://aiida-vasp.readthedocs.org/en/latest>`_
* `Changelog <https://github.com/aiida-vasp/aiida-vasp/blob/develop/CHANGELOG.md>`_

Install and usage:
------------------

Install AiiDA
~~~~~~~~~~~~~

* Install and make sure the `AiiDA prerequisites  <https://aiida.readthedocs.io/projects/aiida-core/en/latest/install/prerequisites.html>`_ work.

* Setup the virtual Python environment as described in the `AiiDA installation instructions <https://aiida.readthedocs.io/projects/aiida-core/en/latest/install/installation.html>`_ and enable it.

* Download latest `AiiDA from GitHub  <https://github.com/aiidateam/aiida_core>`_ by cloning the repository. We suggest to put this into the root folder of the virtual environment. This will be the `aiida_core`.

* Make sure you are on the `develop` branch (should be the default) by issuing in the cloned directory of `aiida_core`::

    git checkout develop

* Install latest `aiida_core` by following the `AiiDA installation instructions <https://aiida.readthedocs.io/projects/aiida-core/en/latest/install/installation.html>`_. In the instruction, the optional packages are descibed. To summarize issue the following in the cloned `aiida_core` directory::

    pip install -e .[option1,option2]

    
Install stable version of the plugin
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A stable version of this plugin is no longer maintained. The current branch `migration_beta` is soon to be merged into `master` and will become our future stable version. Please see next step.

  
Install development version of the plugin
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Download the latest VASP plugin for AiiDA by cloning the 'aiida-vasp <https://github.com/aiida-vasp/aiida-vasp>`_ repository. Again we suggest to put this in the root folder of the virtual environment. This will from here on be referred to as `aiida-vasp`.

* Make sure you are on the `migration_beta` branch (is not the default) by issuing in the cloned directory of `aiida-vasp`::

    git checkout migration_beta

* Install lastest `aiida-vasp` by issuing the following in the cloned `aiida-vasp` directory::

    pip install -e .

* If you want to take part in development and for instance submit a pull request you should issue::

    pip install -e .[dev]


Troubleshoot
~~~~~~~~~~~~

To test wether the installation was successful issue::

   $ verdi plugin aiida.calculation list

   # example output:

   * artithmetic.add
   * templatereplacer
   * vasp.vasp
   * vasp.vasp2w90

You should see vasp.* in the list. There could of course be more plugins listed then what is shown above, depending on your AiiDA plugin stack.

Configure AiiDA
~~~~~~~~~~~~~~~

* See the documentation on how to install `AiiDA plugins <https://aiida.readthedocs.io/projects/aiida-core/en/latest/get_started/index.html>`_ and make sure that you issue the following command, and restart the daemon after `aiida-vasp` have been installed with `pip`::

  reentry scan -r aiida
  verdi daemon restart

* Set up the computers by following the `AiiDA installation guide for computers <https://aiida.readthedocs.io/projects/aiida-core/en/latest/get_started/computers.html>`_. This is typically the login node of you cluster computers. It is also possible to configure a local computer, if you for instance want to run VASP locally (then of course, VASP need to be installed locally).

* Set up the VASP code by following the `AiiDA installation guide for codes <https://aiida.readthedocs.io/projects/aiida-core/en/latest/get_started/codes.html>`_. When you are asked for the plugin, please enter `vasp.vasp`, which is the entry point for the VASP plugin. This is how AiiDA loads plugin related properties. Each code is associated with a computer, so you need to add a separate VASP code entry to every computer.

Uploading potentials
~~~~~~~~~~~~~~~~~~~~

Before using the plugin and VASP, we need to give the plugin access to the potentials. The simples approach to uppload a directory of potentials, say the PAW PBE GW ready potentials you enter the directory where those are present at the local computer. In this directory you would then have the elemental folders, like `Si`, `Si_GW` etc. Then you issue::

  verdi data vasp-potcar uploadfamily --name=some_short_name --description="some_description"

Where `some_short_name` is a name, or a tag that you give this potential family. Typically, it could be `PBE_GW_54` or something similar, depending on the potential you upload. Keep the name rather short, but open for the flexibility to extend with additional families. There in an extra description entry, where `some_description` can be a longer description of this family. By issuing this command the plugin uploads the potentials to the database. The potentials are also protected in VASP by license, so in order to make sure the potentials are not stored for each calculation in the database, only a link to the used potentials is stored, in addition to its SHA512 checksum. The potentials are stored in the database, but not per calculations. If one exports a calculation, there is no potential attached and thus the calculation can be shared with other non-VASP users as well.
  
Using the plugin
~~~~~~~~~~~~~~~~

The plugin is supplied with examples in the `examples` folder. These execute `AiiDA workchains <https://aiida.readthedocs.io/projects/aiida-core/en/latest/working/workflows.html#working-workchains>`_ for a few simple operations that one typically perform using VASP.

The most simple example is the `run_vasp_lean.py` script. This sets up a simple Si cell, with a 9x9x9 k-point grid and some standard INCAR parameters.

To execute this calculation, you need to have configure a computer and a code associated with this computer. Please make sure you have testet the connection to the computer. Then issue::

  python run_vasp_lean.py codename computername --potential-family=some_short_name

The plugin then prepares a VASP calculation and submits its setup to the AiiDA daemon, which again sets this up on the specified `computername` with the core `codename`.

You can then inspect the output by issuing the command::

  verdi process list -a

More details will follow...

