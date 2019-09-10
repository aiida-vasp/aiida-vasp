.. _using_workflows:

Using bundled workflows
=======================
The bundled workflows in AiiDA-VASP is currently only composed of the convergence and bands workflows. The former basically uses the :ref:`converge_workchain` as an entry point to perform the necessary convergence steps for the k-point grid and plane wave cutoff, any relaxations the user need, followed by one final calculation that ejects the parameters of interest. The latter performs first a call to the convergence workflow, followed by a calculation of the standardized special points in reciprocal space. These are combined into one final calculation enabling the extraction of the electronic band structure along these special points in a rather automatic manner. The entry point for the latter is :ref:`master_workchain`. Of course, the convergence workflow can also be entered from the :ref:`master_workchain`.

We would like to encourage users to build workchains and/or compose existing ones into more advanced workflows that we can all share and benefit from.

In the following we will given an example of using the convergence workflow by entering from :ref:`master_workchain`. This would typically illustrate how to perform a `VASP`_ calculation in AiiDA-VASP for a new system, where little is known.

Performing a VASP calculation using the convergence workflow
------------------------------------------------------------

.. _VASP: https://www.vasp.at


