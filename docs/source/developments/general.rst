.. _general_info:

General information
===================

.. note::
   We are currently looking for additional developers. If you are interested, please open an issue on our repository on Github.

Contributing to the development of `AiiDA-VASP`_ is highly encouraged. By contributing you will not only develop a new and more robust way of obtaining your result, but you will also add a very appreciated contribution to the community. The development is a community effort. You are encouraged to use what the community has so far developed. And we would greatly appreciate any contribution or feedback. Please use our issue board and submit a pull request if you have suggestions for code.

If you want to contribute, please have a look at the `open issues`_, notify that you want to contribute on this issue and we will get in touch. If there is no issue for what you want to do, please first open an issue such that we can try to coordinate the effort and make sure we do not double up on the same work.

When developing calculation or workchain plug-ins it should be kept in mind, that to be successfully shared with other researchers, they need to have access to the plug-in as well. This means that any plug-in should ideally be contributed to the official `AiiDA-VASP`_ repository.

Another consideration is that changing a calculation or workchain that has already been used may break provenance, so proceed with extreme caution.

In order to adapt and / or extend the plug-in's calculations and workchains, one should be familiar with the plug-in developer section of the `AiiDA documentation`_.

Also we would encourage users to focus on writing workchains that uses :ref:`converge_workchain` or at least, the very minimal
:ref:`vasp_workchain` as the workchain that is used to fetch `VASP`_ results and as an entry point to the VASP calculation.

.. _AiiDA documentation: https://aiida.readthedocs.io/projects/aiida-core/en/latest/index.html
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _open issues: https://github.com/aiida-vasp/aiida-vasp/issues
.. _VASP: https://www.vasp.at
