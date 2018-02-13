##########
AiiDA-VASP
##########

************
Description
************
`VASP`_ (Vienna Ab initio Simulation Package) is a computer program for 
atomic scale materials modelling, e.g. electronic structure calculations 
and quantum-mechanical molecular dynamics, from first principles.
For more info and a list of features look `here`_. For detailed documentaion
about using VASP take a look at the `documentation page`_ or the `wiki`_

.. _VASP: https://www.vasp.at/
.. _here: https://www.vasp.at/index.php/about-vasp/59-about-vasp
.. _documentation page: https://www.vasp.at/index.php/documentation
.. _wiki: http://cms.mpi.univie.ac.at/wiki/index.php/The_VASP_Manual

*******
Plugins
*******

.. toctree::
   :maxdepth: 3

   VASP <calcs/index>
   Calculatins Base Class <calcs/base>

*********
Workflows
*********

.. toctree::
   :maxdepth: 4

   Autowindows <wf/autowindows>
   Windows <wf/windows>
   Single-Calc <wf/singlecalc>
   Helper <wf/helper>

*****
Usage
*****

Quickstart
==========

For maximum ease of use the python program 'runwf.py' is provided in the aiida_vasp plugin
repository. 

::

   $ python aiida_vasp/utils/runwf.py -h

   usage: runwf.py [-h] [--store-template] workflow input_file

   Run a vasp workflow reading parameters from a json file

   positional arguments:
   workflow          a valid string to load a workflow using WorkflowFactory
   input_file        a .json file containing all necessary parameters for the
                     workflow

   optional arguments:
   -h, --help        show this help message and exit
   --store-template  store an input template in <input file> instead of running
                     the workflow

The general intended way to use it is to first run::

   $ python aiida_vasp/utils/runwf.py --store-template vasp.windows example_input.json
      # stores a json file with input keys and explanations
   $ edit example_input.json
      # edit the json file in your favourite editor, replace all the values with your input
      # parameters
   $ python aiida_vasp/utils/runwf.py vasp.windows example_input.json
      # runs the workflow

Where the vasp.windows commandline parameter would be replaced with the actual
workflow you want to run. If you for example wrote your own workflow and wrote it
analog to the ones included in the plugin, you would use something like::

   $ python aiida_vasp/utils/runwf.py user.<username>.<workflowname>

Assuming you put your workflow into aiida/workflows/user/username/workflowname.py

The runwf.py program is independent of the aiida distribution, so it is possible
to copy and customize it, for example to add automatically
generated categorizing information etc to your wokflows according to your
system of organizing experiments.

It can also be used as a reference on how to read json files into parameter dictionaries
and use those to run workflows.

**********************
Adapting and Extending
**********************

When developping calculation plugins it should be kept in mind, that
for calculations run with them to be successfully shared with other researchers,
those other researchers will need to have access to the calculation plugin as well.
This means that any calculation should be contributed to the aiida_vasp plugin repository.

Another consideration is that changing a calculation that has already been used may
break consistency and reproducability, so proceed with extreme caution.

The same applies in a lesser degree to workflows, as they are not exported from the database
together with the calculation and data nodes. However, care must be taken not to mislead users
by changing the behaviour of a public workflow.

In order to adapt and / or extend the plugin's workflows and calculations, one should be 
familiar with the following parts of the official aiida documentation:

* :doc:`Developping Calculations <aiida:developer_guide/devel_tutorial/code_plugin_int_sum>`
* :doc:`Developping Data Nodes <aiida:developer_guide/devel_tutorial/code_plugin_float_sum>`
* :doc:`Developping Workflows <aiida:developer_guide/workflows>`

The AiiDA-VASP plugin provides some base and helper classes for convenience and to aid with rapid development and / or testing.

* :py:class:`CalcMeta <aiida_vasp.calcs.base.CalcMeta>`, a metaclass that allows writing calculations in a slightly less verbose yet more clear way.
* :py:class:`VaspCalcBase <aiida_vasp.calcs.base.VaspCalcBase>`, takes over common tasks for vasp calculations while providing flexibility using hooks for subclasses.
* :py:class:`BasicCalculation <aiida_vasp.calcs.base.BasicCalculation>`, with functions for writing the basic, common input files.
* :py:class:`BaseParser <aiida.orm.parsers.plugins.vasp.base.BaseParser>`
* :py:class:`VaspParser <aiida.orm.parsers.plugins.vasp.vasp.VaspParser>`
* :py:class:`WorkflowHelper <aiida.workflows.vasp.helper.WorkflowHelper>`, with common convenienct functionality for workflows

There is no base class for workflows, because aiida only recognizes direct subclasses of
:py:class:`Workflow <aiida.workflow.Workflow>` and indirect descendants will raise an
error instead of running.

Calculations
============

For demonstrating how easy small changes are, here is how to implement a Calculation plugin, which takes a ParameterData node for kpoints input to run VASP calculations with path type KPOINTS files.
The first step is to write a new calculation, deriving from BasicCalculation::
   
   from aiida_vasp.calcs.base import BasicCalculation, Input

   class KppathCalculation(BasicCalculation):
      kpoints = Input(types='parameter')

      def write_kpoints(inputdict, destination):
         with open(destination, 'w') as kpoints_file:
            kpoints_file.write(self.get_path_kpoints(self.inp.kpoints.get_dict()))

      def get_path_kpoints(path_spec_dict):
         # KPOINTS_file_content = ...
         # create a string that matches VASP's format
         # for kpoints paths
         return KPOINTS_file_content 

      def check_kpoints(self, kpoints):
         path_spec = kpoints.get_dict()
         # check that the path_spec is in order,
         # raise error otherwise

That's all. 
The line ::

   kpoints = Input(types='parameter')

is automatically converted by :py:class:`CalcMeta <aiida_vasp.calcs.base.CalcMeta>`
into an entry in the classes _use_methods classproperty, overriding the one of BasicCalculation.

The classes write_kpoints method ::

   def write_kpoints(inputdict, destination):
      ...

Overrides the one of BasicCalculation and thus, the new parameter type is handled correctly.

Obviously there would be some work involved in defining the format of the input
ParameterData node and converting it's contents to an appropriate KPOINTS file.

However, there are two problems with this. First of all, it does not take chargedensities or wavefunctions as inputs and second, BasicCalculation's default parser does nothing except verify that an OUTCAR file was created. both can be fixed by simply changing the class so it derives from NscfCalculation::

   ...
   from aiida_vasp.calcs.nscf import NscfCalculation

   class KppathCalculation(NscfCalculation):
   ...

However, it may be educational to see what would be needed to set the extra inputs and default parser by hand::

   class KppathCalculation(BasicCalculation):
      kpoints = Input(types='parameter')
      charge_density = Input(types='vasp.chargedensity')
      wavefunctions = Input(types='vasp.wavefun')
      default_parser = 'vasp.nscf'

      def _prepare_for_submission(self, tempfolder, inputdict):
         calcinfo = super(KppathCalculation, self)._prepare_for_submission(
            tempfoler, inputdict)
         calcinfo.retrieve_list += ['EIGENVAL', 'DOSCAR', 'PROJPAR']

      def write_additional(self, tempfolder, inputdict):
         super(KppathCalculation, self).write_additional(tempfolder, inputdict)
         chgcar_dest = tempfolder.get_abs_path('CHGCAR')
         wavecar_dest = tempfolder.get_abs_path('WAVECAR')
         self.write_chgcar(inputdict, chgcar_dest)
         self.write_wavecar(inputdict, wavecar_dest)

      def write_chgcar(inputdict, destination):
         chgd = inputdict['chargedensity']
         chgcar_path = chgd.get_file_abs_path()
         # copy file at chgcar_path to destination

      def write_wavecar(self, inputdict, destination):
         # analog to write_chgcar

      def write_kpoints(inputdict, destination):
         ... # as before

      ... # as before

      def _init_internal_params(self): ## needed for the default_parser setting to take effect
         super(KppathCalculation, self)._init_internal_params()
         self._update_interal_params()

Here we added more inputs, used write_additional from the superclass to write input files
corresponding to the new input nodes, and extended the _prepare_for_submission from the
superclass (which calls the write_* methods, by the way) to retrieve EIGENVAL and DOSCAR files
from the finished VASP run. The additional output files are then parsed by 'vasp.nscf' into output nodes.

Parsers
-------

Changing a calculation's retrieve_list is only half the way to get additional output nodes,
if none of the current parsers can deal with the additionally retrieved files.

Here's a parser class that delivers the behaviour of the nscf parser used in the above example
but also deals with the PROJCAR file::

   from aiida_vasp.parsers.nscf import NscfParser
   from aiida.orm import DataFactory

   class KppathParser(NscfParser):
      def parse_with_retrieved(self, retrieved):
         super(KppathParser, self).parse_with_retrieved(retrieved)

         proj_file = self.get_file('PROJCAR')
         proj_node = DataFactory('singlefile')(file=proj_file)
         self.add_node('projcar', proj_node)

         return self.result(success=True)

This class makes use of some :py:class:`BaseParser <aiida_vasp.parsers.base.BaseParser>`
convenience functions:

* the superclass keeps a record of all the parsed output nodes added by add_node
* therefore one can just call super.parse_with_retrieved() and in the end return self.result()
  and all the nodes parsed by classes in between are returned as well.
* the get_file call is a convenience function that returns a path to the downloaded (retrieved) file

Naturally in most cases one would not simply like to wrap the output file into a SinglefileData node, but create an array from it, archive it in compressed form, or put it into a datastructure specifically developped for it, however that is documented in :doc:`Developping Data Nodes <aiida:developer_guide/devel_tutorial/code_plugin_float_sum>`.

Workflows
---------

Writing workflows is quite straightforward and explained in detail in :doc:`Developping Workflows <aiida:developer_guide/workflows>`.
There is no base class provided, since all workflows must be derived directly from aiida's Workflow class.

However there is :py:class:`WorkflowHelper <aiida.workflows.vasp.helper.WorkflowHelper>`.
The provided workflows make use of it as a common way to provide the possibility to
document the parameters they use, especially the common ones. It also allows a unified
way to verify parameter consistency as well as some common logging functionality.
