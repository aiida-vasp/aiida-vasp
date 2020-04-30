.. _run_vasp_builder:

=======================================
5. Running VASP calculations explicitly
=======================================

We will now try to build a simple call script that executes a workchain,
basically a of one-shot VASP calculation for wurtzite-AlN. When that is done,
we will modify this example is to also use the relax
workchain, which enables relaxations of the structure.


Before running VASP calculation
--------------------------------

It is recommended to well symmetrize the crystal strucutre of interest
if the space group type is known. This can be done by using
spglib. Spglib with python interface can be installed either pip or
conda, e.g.::

   % pip install spglib

or::

   % conda install -c conda-forge spglib

The usage of spglib is found at
https://atztogo.github.io/spglib/. Let's take an example of POSCAR of
wurtzite-type SiC structure, which can be obtained from the Materials
project database, https://materialsproject.org/materials/mp-7140/

::

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

This is already symmetrized but we can symmetrize more if we want more
number of digits for :math:`\sqrt{3}` and :math:`1/3` that often
appear in hexagonal crystals. An ad-hoc script to symmetrize this may
be like

::

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

The space group type is found by
``print(spglib.get_spacegroup(cell))``, which should give ``P6_3mc
(186)``  in this example, and we have to be sure at least that this is
the correct one. Then, with::

   np.set_printoptions(precision=15)
   [print(sym_cell[i]) for i in range(3)]

we see::

   [[ 3.092007293580808  0.                 0.               ]
    [-1.546003646790404  2.677756864927749  0.               ]
    [ 0.                 0.                 5.073347         ]]
   [[0.333333333333333 0.666666666666667 0.499589         ]
    [0.666666666666667 0.333333333333333 0.999589         ]
    [0.333333333333333 0.666666666666667 0.875411         ]
    [0.666666666666667 0.333333333333333 0.375411         ]]
   [14 14  6  6]

The third value of the atomic position, e.g., ``0.499589``, may be
expected to be zero or 0.5 by feeling, but this is not possible to be
done by spglib, because it has freedome against rigid shift along c
that doesn't alter the symmetry.


Use of AiiDA-VASP
-----------------

A typical script to launch a VASP caluculation (wurtzite-type SiC) is
something like::

   import numpy as np
   from aiida.common.extendeddicts import AttributeDict
   from aiida.manage.configuration import load_profile
   from aiida.orm import Bool, Str, Code, Int, Float
   from aiida.plugins import DataFactory, WorkflowFactory
   from aiida.engine import submit
   load_profile()


   def launch_aiida(structure, code_string, resources,
                    label="SiC VASP calculation"):
       Dict = DataFactory('dict')
       KpointsData = DataFactory("array.kpoints")

       incar_dict = {
           'PREC': 'Accurate',
           'IBRION': -1,
           'EDIFF': 1e-8,
           'NELMIN': 5,
           'NELM': 100,
           'ENCUT': 500,
           'IALGO': 38,
           'ISMEAR': 0,
           'SIGMA': 0.01,
           'GGA': 'PS',
           'LREAL': False,
           'LCHARG': False,
           'LWAVE': False,
       }

       kpoints = KpointsData()
       kpoints.set_kpoints_mesh([6, 6, 4], offset=[0, 0, 0.5])

       options = {'resources': resources,
                  'account': 'nn9995k',
		  'max_memory_kb': 1024000,
                  'max_wallclock_seconds': 3600 * 10}

       potential_family = 'pbe'
       potential_mapping = {'Si': 'Si', 'C': 'C'}

       parser_settings = {'add_energies': True,
                          'add_forces': True,
                          'add_stress': True}

       code = Code.get_from_string(code_string)
       Workflow = WorkflowFactory('vasp.vasp')
       builder = Workflow.get_builder()
       builder.code = code
       builder.parameters = Dict(dict=incar_dict)
       builder.structure = structure
       builder.settings = Dict(dict={'parser_settings': parser_settings})
       builder.potential_family = Str(potential_family)
       builder.potential_mapping = Dict(dict=potential_mapping)
       builder.kpoints = kpoints
       builder.options = Dict(dict=options)
       builder.metadata.label = label
       builder.metadata.description = label
       builder.clean_workdir = Bool(False)

       node = submit(builder)
       return node


   def get_structure_SiC():
       """Set up SiC cell

       Si C
          1.0
            3.0920072935808083    0.0000000000000000    0.0000000000000000
           -1.5460036467904041    2.6777568649277486    0.0000000000000000
            0.0000000000000000    0.0000000000000000    5.0733470000000001
        Si C
          2   2
       Direct
          0.3333333333333333  0.6666666666666665  0.4995889999999998
          0.6666666666666667  0.3333333333333333  0.9995889999999998
          0.3333333333333333  0.6666666666666665  0.8754109999999998
          0.6666666666666667  0.3333333333333333  0.3754109999999997

       """

       StructureData = DataFactory('structure')
       a = 3.092
       c = 5.073
       lattice = [[a, 0, 0],
                  [-a / 2, a / 2 * np.sqrt(3), 0],
                  [0, 0, c]]
       structure = StructureData(cell=lattice)
       for pos_direct, symbol in zip(
               ([1. / 3, 2. / 3, 0],
                [2. / 3, 1. / 3, 0.5],
                [1. / 3, 2. / 3, 0.375822],
                [2. / 3, 1. / 3, 0.875822]), ('Si', 'Si', 'C', 'C')):
           pos_cartesian = np.dot(pos_direct, lattice)
           structure.append_atom(position=pos_cartesian, symbols=symbol)
       return structure


   def main(code_string, resources):
       structure = get_structure_SiC()
       launch_aiida(structure, code_string, resources)


   if __name__ == '__main__':
       code_string = 'vasp@saga'
       resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 20}
       main(code_string, resources)

Once the calculation is done, we locate the ``<pk>`` by using::

  verdi process list

Pick the most recent ``VaspWorkChain`` process and then we can watch the results using
for instance the verdi shell::

  verdi shell

And then we load the node::

   In [1]: n = load_node(<pk>)

   In [2]: n.outputs.energies.get_array('energy_no_entropy')
   Out[2]: array([-31.80518222])

   In [3]: n.outputs.stress.get_array('final')
   Out[3]:
   array([[-29.89502712,   0.        ,   0.        ],
          [  0.        , -29.89502712,   0.        ],
          [  0.        ,   0.        , -29.47075517]])

When we want to fully relax a crystal structure, the above script is
modified as follows:

1. Replace ``WorkflowFactory('vasp.vasp')`` by ``WorkflowFactory('vasp.relax')``
2. Remove ``IBRION`` from ``incar_dict``
3. Add the following setting::

       relax = AttributeDict()
       relax.perform = Bool(True)        # Turn on relaxation of the structure
       relax.force_cutoff = Float(1e-5)  # Relax force cutoff
       relax.steps = Int(10)             # Relax number of ionic steps
       relax.positions = Bool(True)      # Relax atomic positions
       relax.shape = Bool(True)          # Relax cell shape (alpha, beta, gamma)
       relax.volume = Bool(True)         # Relax volume
       builder.relax = relax
       builder.verbose = Bool(True)

The lattice parameters of the relax crystal structure is found by

::

   In [1]: n = load_node(<PK>)

   In [2]: n.outputs.structure_relaxed.cell
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
