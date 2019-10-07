.. _developments:

======================
Adapting and Extending
======================

When developing calculation plug-ins it should be kept in mind, that
for calculations run with them to be successfully shared with other researchers,
those other researchers will need to have access to the calculation plug-in as well.
This means that any calculation should be contributed to the ``github.com/aiida-vasp/aiida-vasp`` plug-in repository.

Another consideration is that changing a calculation that has already been used may
break consistency and reproducibility, so proceed with extreme caution.

The same applies in a lesser degree to workflows, as they are not exported from the database
together with the calculation and data nodes. However, care must be taken not to mislead users
by changing the behaviour of a public workflow.

In order to adapt and / or extend the plug-in's workflows and calculations, one should be
familiar with the following parts of the official AiiDA documentation:

* :doc:`Developping Calculations <aiida:developer_guide/devel_tutorial/code_plugin_int_sum>`
* :doc:`Developping Data Nodes <aiida:developer_guide/devel_tutorial/code_plugin_float_sum>`
* :doc:`Developping Workflows <aiida:developer_guide/workflows>`

The AiiDA-VASP plug-in provides some base and helper classes for convenience and to aid with rapid development and / or testing.

* :py:class:`CalcMeta <aiida_vasp.calcs.base.CalcMeta>`, a metaclass that allows writing calculations in a slightly less verbose yet more clear way.
* :py:class:`VaspCalcBase <aiida_vasp.calcs.base.VaspCalcBase>`, takes over common tasks for VASP calculations while providing flexibility using hooks for subclasses.
* :py:class:`BasicCalculation <aiida_vasp.calcs.base.BasicCalculation>`, with functions for writing the basic, common input files.
* :py:class:`BaseParser <aiida.orm.parsers.plugins.vasp.base.BaseParser>`
* :py:class:`VaspParser <aiida.orm.parsers.plugins.vasp.vasp.VaspParser>`
