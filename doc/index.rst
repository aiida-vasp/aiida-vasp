VASP
----

Description
^^^^^^^^^^^
`VASP`_ (Vienna Ab initio Simulation Package) is a computer program for 
atomic scale materials modelling, e.g. electronic structure calculations 
and quantum-mechanical molecular dynamics, from first principles.
For more info and a list of features look `here`_. For detailed documentaion
about using VASP take a look at the `documentation page`_ or the `wiki`_

.. _VASP: https://www.vasp.at/
.. _here: https://www.vasp.at/index.php/about-vasp/59-about-vasp
.. _documentation page: https://www.vasp.at/index.php/documentation
.. _wiki: http://cms.mpi.univie.ac.at/wiki/index.php/The_VASP_Manual

Plugins
^^^^^^^

General Purpose
===============
These are calculation plugins which aim to expose the full flexibility
of the associated programs and will retrieve all output files, making
no assumption about which output will be needed.
This comes at the cost that potentially big files that won't be used
in analysis will be retrieved and permanently stored in your data base.
Therefore it is recommended to use these only where no special purpose
calculation is available for the task.

.. toctree::
   :maxdepth: 4

   Vasp5Caclculation <vasp5>
   WannierCalculation <wannier>

Special Purpose
===============                            
These calculations are used for specific use cases, like generating
CHGCAR and WAVECAR files or band structures, DOS, or getting input files
for wannier90.x calculations.

They are mostly designed to be used in conjunction with other calculations
in a workflow.

.. toctree::
   :maxdepth: 4

   ScfCalculation - Selfconsisten Run <scf>
   NscfCalculcation - Nonselfconsisten Run <nscf>
   AmnCalculation - Vasp2Wannier90 Projections <amn>
