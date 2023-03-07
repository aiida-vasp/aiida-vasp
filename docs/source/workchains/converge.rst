.. _converge_workchain:

==================
Converge workchain
==================

The ``ConvergeWorkChain`` is intended to be used to determine the appropriate energy cutoff and reciprocal space mesh as to converge the ground state total energy. One can control whether only the energy, reciprocal mesh or both parameters are used.

When performing the energy cutoff convergence if no k-point mesh is given the energy cutoff convergence tests will be done with an auto-generated k-point mesh. The distance between the points is either user provided or default values are used. Similarly, for the reciprocal mesh, if a provided energy cutoff has been given that will be used, otherwise if a cutoff value has been found from the energy convergence that would be used instead.

.. warning::
   Please not that this documentation is not up to date. Please consult the example ``run_converge.py`` and the
   documentation in the code of ``ConvergeWorkChain``.
