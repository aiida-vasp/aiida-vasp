=================
AiiDA Vasp Plugin
=================

.. image:: https://travis-ci.org/aiida-team/aiida-vasp.svg?branch=develop
    :target: https://travis-ci.org/aiida-team/aiida-vasp

.. image:: https://readthedocs.org/projects/aiida-vasp/badge/?version=latest
   :target: http://aiida-vasp.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://coveralls.io/repos/github/aiida-team/aiida-vasp/badge.svg?branch=coveralls
   :target: https://coveralls.io/github/aiida-team/aiida-vasp?branch=coveralls

This is a plugin to `AiiDA <www.aiida.net/?page_id=264>`_ to run calculations with the ab-initio program `VASP <https://www.vasp.at/>`_.

* `Documentation <https://aiida-vasp.readthedocs.org/en/latest>`_
* `Changelog <https://github.com/aiida-team/aiida-vasp/blob/develop/CHANGELOG.md>`_

Install and usage:
------------------

Install AiiDA
~~~~~~~~~~~~~

* Download, install and setup -> `aiida's documentation <http://aiida-core.readthedocs.org/en/stable/>`_

Install this plugin
~~~~~~~~~~~~~~~~~~~

From the aiida-vasp folder use::

   $ pip install .
   # or
   $ pip install -e . # <-- for plugin developers

To test wether the installation was successful use::

   $ verdi calculation plugins 

   # example output:

   ## Pass as a further parameter one (or more) plugin names
   ## to get more details on a given plugin.
   * codtools.cifcellcontents
   * codtools.cifcodcheck
   * codtools.cifcodnumbers
   * codtools.ciffilter
   * codtools.cifsplitprimitive
   * quantumespresso.cp
   * quantumespresso.pw
   * quantumespresso.pwimmigrant
   * simpleplugins.templatereplacer
   * vasp.amn
   * vasp.scf
   * vasp.nscf
   * vasp.vasp
   * vasp.wannier

You should see vasp.* in the list

Configure the code
~~~~~~~~~~~~~~~~~~

See `aiida docs <http://aiida-core.readthedocs.org/en/stable/setup/computerandcodes.html#computer-setup-and-configuration>`_
on how to set up computers and codes. Note that you need at least one computer configured and a VASP executable on it
in order to use this plugin.

Using this plugin
~~~~~~~~~~~~~~~~~

For a tryout session, use the script ``run_vasp.py`` in the examples folder. Assuming you have already set up a code and computer, usage is as follows::

   $ python run_vasp.py --help  # shows help message
   $ # first time you need to import the potentials
   $ python run_vasp.py --import-from=<path-to-potentials-folder> --queue=<your-compute-queue> <code> <computer>
     ... will take a while the first time due to importing potentials
   $ # from  now on you can simply run
   $ # python run_vasp.py --queue=<queue> <code> <computer>
   $ verdi calculation list
     .. should contain a vasp calculation


