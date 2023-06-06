.. _workflows:

Workflows
=========
By composing one or several `Workchain`_ classes, one can make a workflow. A workflow can of course only consits of one rather trivial operation, but is becoming truly usefull when the task at hand is complicated and time consuming to construct, inspect, nourish and analyze by hand. With the introduction of high-throughput calculations and more composed single calculations using `VASP`_ such a feature is welcoming.

The bundled workflows in AiiDA-VASP is currently only composed of the convergence and bands workflow. The former basically uses the :ref:`converge_workchain` as an entry point to perform the necessary convergence steps for the k-point grid and plane wave cutoff, any relaxations the user need, followed by one final calculation that ejects the parameters of interest. The latter performs first a call to the convergence workflow, followed by a calculation of the standardized special points in reciprocal space. These are combined into one final calculation enabling the extraction of the electronic band structure along these special points in a rather automatic manner.

We would like to encourage users to build workchains and/or compose existing ones into more advanced workflows that we can all share and benefit from.

.. _VASP: https://www.vasp.at
.. _Workchain: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/workflows.html#work-chains
