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

