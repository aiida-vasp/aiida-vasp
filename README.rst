=================
AiiDA Vasp Plugin
=================

Install and usage:
------------------

Install AiiDA
~~~~~~~~~~~~~

1. Download from [www.aiida.net/?page_id=264](www.aiida.net -> Download) (only EPFL version supports ssh)
2. install and setup -> [http://aiida-core.readthedocs.org/en/stable/](aiida's documentation)

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

See [http://aiida-core.readthedocs.org/en/stable/setup/computerandcodes.html#computer-setup-and-configuration](aiida docs)
on how to set up computers and codes. Note that you need at least one computer configured and a VASP executable on it
in order to use this plugin.

Using this plugin
~~~~~~~~~~~~~~~~~

This section is still in construction
