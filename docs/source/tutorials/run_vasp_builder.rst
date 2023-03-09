.. _run_vasp_builder:

=======================
1. Structure relaxation
=======================

In this tutorial, we will first perform a calculation of the electronic
structure of wurtzite-AlN, following by a relaxation of the structure.
The idea here is to understand how we can use the ``builder`` to set up
the inputs for the call script and how to interface with the bundled
workchain for relaxing structures.

Before continuing it is worthwile to consider symmetry and work with that in a more
general manner. It is good practice to familiar ourself about these topics from
time to time as much of it is easily forgotten and taken for granted, giving us
problems downstream that could be avoided with a bit better oversight initially.

Symmetrizing structures
-----------------------

It is recommended to symmetrize the crystal structure of interest
if the space group type is known. This can be done by using for instance
`spglib`_, which should already be installed when installing `AiiDA-VASP`_.
Let's take an example ``POSCAR`` of wurtzite-type SiC structure,
which can be `obtained`_ from the Materials project database::

  Si2 C2
  1.0
  3.092007 0.000000 0.000000
  -1.546004 2.677757 0.000000
  0.000000 0.000000 5.073347
  Si C
  2 2
  direct
  0.333333 0.666667 0.499589 Si
  0.666667 0.333333 0.999589 Si
  0.333333 0.666667 0.875411 C
  0.666667 0.333333 0.375411 C

This is already symmetrized but we can symmetrize it even more if we want more
number of digits for :math:`\sqrt{3}` and :math:`1/3` that often
appear in hexagonal crystals. An ad-hoc Python script to symmetrize this may look like::

  import numpy as np
  import spglib

  poscar_lines = """Si2 C2
  1.0
  3.092007 0.000000 0.000000
  -1.546004 2.677757 0.000000
  0.000000 0.000000 5.073347
  Si C
  2 2
  direct
  0.333333 0.666667 0.499589 Si
  0.666667 0.333333 0.999589 Si
  0.333333 0.666667 0.875411 C
  0.666667 0.333333 0.375411 C""".splitlines()

  lattice = np.genfromtxt(poscar_lines[2:5]).reshape(3, 3)
  points = np.genfromtxt(poscar_lines[8:12]).reshape(4, -1)[:, :3]
  numbers = [14, 14, 6, 6]
  cell = (lattice, points, numbers)

  sym_cell = spglib.refine_cell(cell)
  print('Spacegroup:', spglib.get_spacegroup(cell))
  np.set_printoptions(precision=15)
  [print(sym_cell[i]) for i in range(3)]

Upon saving and executing this Python script, the space group should be printed as ``P6_3mc (186)``.
And the symmetrized cell should be::

   [[ 3.092007293580808  0.                 0.               ]
    [-1.546003646790404  2.677756864927749  0.               ]
    [ 0.                 0.                 5.073347         ]]
   [[0.333333333333333 0.666666666666667 0.499589         ]
    [0.666666666666667 0.333333333333333 0.999589         ]
    [0.333333333333333 0.666666666666667 0.875411         ]
    [0.666666666666667 0.333333333333333 0.375411         ]]
   [14 14  6  6]

The third component of the atomic position, i.e. ``0.499589``, may be
expected to be zero or 0.5 by feeling, but this can never be achieved exactly by ``spglib``,
or other codes working with numerics. The take home message from this is that
you should know what kind of structure you are looking into and from time to time make checks
to make sure you are still where you think you are. Let us continue and relax this structure
using `VASP`_ and the plugin.

Relaxation of the structure
---------------------------

Here we will relax the structure, but first, we will perform an initial run without relaxing the structure
so that you get familiar with the system and also are able to make the quick changes to enable relaxation from
the calling script.

#. First assemble a script to launch a VASP caluculation (wurtzite-type SiC). A suggestion would be::

     DDEJIDE

   And save it to for instance ``run_vasp_sic.py``.

#. Run the saved script to launch the calculation::

     $ python run_vasp_sic.py

#. Wait a bit and once the ``VaspWorkChain`` is finalized, get its ``PK`` using ``verdi process list -a``, here ``2476``.

#. Check the attached nodes of ``2476``::



#. And inspect for instance the ``energies`` output node::



   Do not take these values for granted and compare them to yours. They depend on the
   system you executed, potential used etc.

#. Start ``verdi shell`` and then load the node::

     In [1]: n = load_node(<PK>)

     In [2]: n.outputs.energies.get_array('energy_extrapolated')
     Out[2]: array([-31.80518222])

     In [3]: n.outputs.stress.get_array('final')
     Out[3]:
     array([[-29.89502712,   0.        ,   0.        ],
     [  0.        , -29.89502712,   0.        ],
     [  0.        ,   0.        , -29.47075517]])

   Let us now modify the script so that we perform a structure relaxation.
   If we want to fully relax the crystal structure, we need to modify the script accordingly.

#. Exit ``verdi shell`` by typing ``exit``.

#. Replace ``WorkflowFactory('vasp.vasp')`` with ``WorkflowFactory('vasp.relax')``

#. Remove the ``IBRION`` entry from ``incar_dict``

#. Add add after ``builder.clean_workdir = Bool(False)`` the following::

     relax = AttributeDict()
     relax.perform = Bool(True)        # Turn on relaxation of the structure
     relax.force_cutoff = Float(1e-5)  # Relax force cutoff
     relax.steps = Int(10)             # Relax number of ionic steps
     relax.positions = Bool(True)      # Relax atomic positions
     relax.shape = Bool(True)          # Relax cell shape (alpha, beta, gamma)
     relax.volume = Bool(True)         # Relax volume
     builder.relax = relax
     builder.verbose = Bool(True)

   The plugin will then set the correct ``ISIF`` and ``EDIFFG`` etc. The point of using dedicated
   settings like this is twofolded: (i) it makes more sense to the user, and (ii) it makes this
   workflow independent on `VASP`_ and can in principle be executed with any other backend, say
   Quantum Espresso as long as the conversion in the backend is done properly. The long term
   goal of the development of these plugins is that we will eventually have a more unified interface
   for the workflows that in principle can be code independent.

#. Save the modified script and relaunch it.

#. Locate the ``PK`` of the finalized ``RelaxWorkChain`` and launch the ``verdi shell`` again.
   Then locate the relaxed structure and the stress::

     In [1]: n = load_node(<PK>)

     In [2]: n.outputs.relax__structure.cell
     Out[2]:
     [[3.07798535, 0.0, 0.0],
     [-1.53899268, 2.66561351, 0.0],
     [0.0, 0.0, 5.04931673]]

     In [3]: n.outputs.stress.get_array('final')
     Out[3]:
     array([[-0.01708304,  0.        ,  0.        ],
     [ 0.        , -0.01708304,  0.        ],
     [ 0.        ,  0.        , -0.00809151]])

There are more options for the relax workchain, e.g., running VASP
several time iteratively until convergence, which is used in the bulk
modulus example in the next section.

After the relaxation, sometimes the crystal symmetry can be slightly
broken by the VASP calculation, especially for hexagonal crystals. So
it is recommended to symmetrize the final structure if this is the case.


.. _obtained: https://materialsproject.org/materials/mp-7140/
.. _spglib: https://spglib.github.io/spglib/
.. _VASP: https://www.vasp.at
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
