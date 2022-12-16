.. _wannier_calculation:

=============================
VASP to Wannier90 calculation
=============================

Description
-----------

Works exactly the same as :ref:`vasp_calculation`, from which it is derived, with the difference
that it takes a ``wannier_parameters`` input as well, for use with the VASP2Wannier90 interface. This requires a VASP code compiled with ``LWANNIER90=True``.

Inputs
------

* :ref:`parameters <vasp-input-parameters>`
* :ref:`kpoints <vasp-input-kpoints>`
* :ref:`structure <vasp-input-structure>`
* :ref:`potential <vasp-input-potential>`
* :ref:`charge_density <vasp-input-charge>`
* :ref:`wavefunctions <vasp-input-wave>`
* :ref:`wannier_parameters <vasp-input-wannier_parameters>`

Outputs
-------

* :ref:`results <vasp-output-results>`
* :ref:`bands <vasp-output-bands>`
* :ref:`dos <vasp-output-dos>`
* :ref:`wannier_parameters <vasp-output-wannier_parameters>`
* :ref:`wannier_data <vasp-output-wannier_data>`
