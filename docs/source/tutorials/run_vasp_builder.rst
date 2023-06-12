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
expected to be zero or 0.5, but this can never be achieved exactly by ``spglib``,
or other codes working with numerics. The take home message from this is that
you should know what kind of structure you are looking into and from time to time make checks
to make sure you are still where you think you are. Let us continue and relax this structure
using `VASP`_ and the plugin.

Relaxation of the structure
---------------------------

Here we will relax the structure, but first, we will perform an initial run without relaxing the structure
so that you get familiar with the system and also are able to make the quick changes to enable relaxation from
the calling script.

Initial static run
^^^^^^^^^^^^^^^^^^

#. First assemble a script to launch a VASP caluculation using the wurtzite-type SiC structure
   . The scripts to launch certain calculations can be designed in many different way.
   Let us fetch an example `AiiDA-VASP`_ run file::

     $ wget https://github.com/aiida-vasp/aiida-vasp/raw/master/tutorials/run_sic.py

#. Inspect the file, which has the following content:

   .. literalinclude:: ../../../tutorials/run_sic.py

#. Change the relevant bits, e.g. the ``code_string``, the ``options`` and the ``potential_family`` for your
   stored setup. Please consult previous tutorials for details on this.

#. Run the modified script to launch the calculation::

     $ python run_sic.py

#. Wait a bit and once the ``VaspWorkChain`` is finalized, get its ``PK`` using ``verdi process list -a``, here ``2663``.

#. Check the output nodes of ``2663``::

     $ verdi process show 2663
     Property     Value
     -----------  ------------------------------------
     type         VaspWorkChain
     state        Finished [0]
     pk           2663
     uuid         cc1f4f29-9c96-4e17-ba2f-3a140f789d49
     label        SiC VASP calculation
     description  SiC VASP calculation
     ctime        2023-03-22 10:47:17.841728+01:00
     mtime        2023-03-22 10:49:45.096366+01:00

     Inputs               PK  Type
     -----------------  ----  -------------
     clean_workdir      2660  Bool
     code                  3  InstalledCode
     kpoints            2658  KpointsData
     max_iterations     2661  Int
     options            2659  Dict
     parameters         2653  Dict
     potential_family   2656  Str
     potential_mapping  2657  Dict
     settings           2655  Dict
     structure          2654  StructureData
     verbose            2662  Bool

     Outputs          PK  Type
     -------------  ----  ----------
     energies       2668  ArrayData
     forces         2670  ArrayData
     misc           2669  Dict
     remote_folder  2666  RemoteData
     retrieved      2667  FolderData
     stress         2671  ArrayData

     Called          PK  Type
     ------------  ----  ---------------
     iteration_01  2665  VaspCalculation

     Log messages
     ---------------------------------------------
     There are 3 log messages for this calculation
     Run 'verdi process report 2663' to see them


#. And inspect for instance the ``energies`` output node::

     $ verdi data array show 2668
     {
	 "electronic_steps": [
	     1
	 ],
	 "energy_extrapolated": [
	     -30.09913368
	 ],
	 "energy_extrapolated_electronic": [
	     -30.09913368
	 ]
     }

   Do not take these values for granted and compare them to yours. They depend on the
   system you executed, potential used etc.

#. We can also inspect it from Python. In order to make all the AiiDA machinery
   available in Python, we can use ``verdi shell``. Start ``verdi shell`` and then load the node::


     $ verdi shell
     noPython 3.10.10 (main, Mar  5 2023, 22:26:53) [GCC 12.2.1 20230201]
     Type 'copyright', 'credits' or 'license' for more information
     IPython 7.34.0 -- An enhanced Interactive Python. Type '?' for help.

     In [1]: n = load_node(2663)

     In [2]: n.outputs.energies.get_array('energy_extrapolated')
     Out[2]: array([-30.09913368])

     In [3]: n.outputs.stress.get_array('final')
     Out[3]:
     array([[-0.01000677,  0.        ,  0.        ],
	    [-0.        , -0.01000677,  0.        ],
	    [ 0.        ,  0.        ,  0.34873446]])

#. Exit ``verdi shell`` by typing ``exit``.

Now that we have verified that the script works, let us extend it to enable relaxation of the
structure.

Relaxation run
^^^^^^^^^^^^^^

Let us now modify the script so that we perform a structure relaxation.
If we want to fully relax the crystal structure, we need to modify the script accordingly.

#. Open the ``run_sic.py`` launch script again.

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

   In the end, the necessary changes to ``run_sic.py`` can be summarize as:

   .. literalinclude:: ../../../tutorials/run_sic_relax.py
      :diff: ../../../tutorials/run_sic.py

   The plugin will then set the correct ``ISIF`` and ``EDIFFG`` etc. The point of using dedicated
   settings like this is twofolded: (i) it makes more sense to the user, and (ii) it makes this
   workflow independent on `VASP`_ and can in principle be executed with any other backend, say
   Quantum Espresso as long as the conversion in the backend is done properly. This is the role of the
   :py:class:`ParametersMassage<aiida_vasp.assistant.parameters.ParametersMassage>`. The long term
   goal of the development of these plugins is that we will eventually have a more unified interface
   for the workflows that in principle can be code independent.

#. Save the modified script and relaunch it.

#. Locate the ``PK`` of the finalized ``RelaxWorkChain``, in this case ``2698``
   and launch ``verdi shell`` again. Then locate the relaxed structure and the stress::


     $ verdi shell
     Python 3.10.10 (main, Mar  5 2023, 22:26:53) [GCC 12.2.1 20230201]
     Type 'copyright', 'credits' or 'license' for more information
     IPython 7.34.0 -- An enhanced Interactive Python. Type '?' for help.

     In [1]: n = load_node(2698)

     In [2]: n.outputs.relax.structure.cell
     Out[2]:
     [[3.09208384, 0.0, 0.0],
      [-1.54604192, 2.67782316, 0.0],
      [0.0, 0.0, 5.07299923]]

     In [3]: n.outputs.stress.get_array('final')
     Out[3]:
     array([[-0.00109338,  0.        ,  0.        ],
	    [ 0.        , -0.00109338,  0.        ],
	    [ 0.        ,  0.        , -0.00031531]])

There are more options for the relax workchain, e.g., running VASP
several time iteratively until convergence, which is used in the bulk
modulus example in the next section.

After the relaxation, sometimes the crystal symmetry can be slightly
broken by the VASP calculation, especially for hexagonal crystals. So
it is recommended to symmetrize the final structure if this is the case, depending
on what you want to use it for.

.. _obtained: https://materialsproject.org/materials/mp-7140/
.. _spglib: https://spglib.github.io/spglib/
.. _VASP: https://www.vasp.at
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
