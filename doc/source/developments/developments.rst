.. _developments:

======================
Adapting and Extending
======================

When developing calculation or workchain plug-ins it should be kept in mind, that to be successfully shared with other researchers,
they need to have access to the plug-in as well.
This means that any plug-in should ideally be contributed to the official `AiiDA-VASP`_ repository.

Another consideration is that changing a calculation or workchain that has already been used may
break provenance, so proceed with extreme caution.

In order to adapt and / or extend the plug-in's calculations and workchains, one should be
familiar with the plug-in developer section of the `AiiDA documentation`_.

Also we would encourage users to focus on writing workchains that uses :ref:`converge_workchain` or at least, the very minimal
:ref:`vasp_workchain` as the workchain that is used to fetch `VASP`_ results.

.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _AiiDA documentation: https://aiida.readthedocs.io/projects/aiida-core/en/latest/index.html
.. _VASP: https://www.vasp.at/
