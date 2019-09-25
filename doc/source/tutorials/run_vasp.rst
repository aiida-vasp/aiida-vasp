.. _run_vasp:

========================
Running VASP calculation
========================

An example of one-shot VASP calculation for wurtzite-AlN is
given. Then this example is slightly modified to use the relax
workchain.


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
wurtzite-type AlN structure, which can be obtained from the Materials
project database, https://materialsproject.org/materials/mp-661/

::

   Al2 N2
   1.0
   3.128588 0.000000 0.000000
   -1.564294 2.709437 0.000000
   0.000000 0.000000 5.016955
   Al N
   2 2
   direct
   0.333333 0.666667 0.999287 Al
   0.666667 0.333333 0.499287 Al
   0.333333 0.666667 0.380713 N
   0.666667 0.333333 0.880713 N

This is already symmetrized but we can symmetrize more if we want more
number of digits for :math:`\sqrt{3}` and :math:`1/3` that often
appear in hexagonal crystals. An ad-hoc script to symmetrize this may
be like

::

   import numpy as np
   import spglib

   poscar_lines = """ Al2 N2
   1.0
   3.128588 0.000000 0.000000
   -1.564294 2.709437 0.000000
   0.000000 0.000000 5.016955
   Al N
   2 2
   direct
   0.333333 0.666667 0.999287 Al
   0.666667 0.333333 0.499287 Al
   0.333333 0.666667 0.380713 N
   0.666667 0.333333 0.880713 N""".splitlines()

   lattice = np.genfromtxt(poscar_lines[2:5]).reshape(3, 3)
   points = np.genfromtxt(poscar_lines[8:12]).reshape(4, -1)[:, :3]
   numbers = [13, 13, 7, 7]
   cell = (lattice, points, numbers)

   sym_cell = spglib.refine_cell(cell)

The space group type is found by
``print(spglib.get_spacegroup(cell))``, which should give ``P6_3mc
(186)``  in this example, and we have to be sure at least that this is
the correct one. Then, with::

   np.set_printoptions(precision=15)
   [print(sym_cell[i]) for i in range(3)]

we see::

   [[ 3.128588135976751  0.                 0.               ]
    [-1.564294067988376  2.70943680373447   0.               ]
    [ 0.                 0.                 5.016955         ]]
   [[0.333333333333333 0.666666666666667 0.999287         ]
    [0.666666666666667 0.333333333333333 0.499287         ]
    [0.333333333333333 0.666666666666667 0.380713         ]
    [0.666666666666667 0.333333333333333 0.880713         ]]
   [13 13  7  7]

The third value of the atomic position, e.g., ``0.999287``, may be
expected to be zero by feeling, but this is not possible to be done by
spglib, because it has freedome against rigid shift along c that
doesn't alter the symmetry.


Use of AiiDA-VASP
-----------------

A typical script to launch a VASP caluculation (rutile-type SnO2) is
something like::

   import numpy as np
   from aiida.manage.configuration import load_profile
   from aiida.orm import Bool, Float, Int, Str, Code
   from aiida.plugins import DataFactory, WorkflowFactory
   from aiida.engine import submit
   load_profile()


   def launch_aiida(structure, code_string, resources,
                    label="SnO2 VASP calculation"):
       Dict = DataFactory('dict')
       KpointsData = DataFactory("array.kpoints")

       incar_dict = {
           'PREC': 'Accurate',
           'IBRION': -1,
           'EDIFF': 1e-5,
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
       kpoints.set_kpoints_mesh([4, 4, 6], offset=[0.5, 0.5, 0.5])

       options = {'resources': resources,
                  'max_wallclock_seconds': 3600 * 10}

       potential_family = 'PBE.54'
       potential_mapping = {'Sn': 'Sn', 'O': 'O'}

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


   def get_structure_SnO2():
       """Set up SnO2 structure

       SnO2
          1.0
            4.77 0.00 0.00
            0.00 4.77 0.00
            0.00 0.00 3.22
        Sn O
          2 4
       Direct
          0.000 0.000 0.000
          0.500 0.500 0.500
          0.306 0.306 0.000
          0.694 0.694 0.000
          0.194 0.806 0.500
          0.806 0.194 0.500

       """

       StructureData = DataFactory('structure')
       a = 4.77
       c = 3.22
       lattice = [[a, 0, 0],
                  [0, a, 0],
                  [0, 0, c]]
       structure = StructureData(cell=lattice)
       u = 0.306
       for pos_direct, symbol in zip(
               ([0, 0, 0],
                [0.5, 0.5, 0.5],
                [u, u, 0],
                [1 - u, 1 - u, 0],
                [0.5 - u, 0.5 + u, 0.5],
                [0.5 + u, 0.5 - u, 0.5]), ('Sn', 'Sn', 'O', 'O', 'O', 'O')):
           pos_cartesian = np.dot(pos_direct, lattice)
           structure.append_atom(position=pos_cartesian, symbols=symbol)
       return structure


   def main(code_string, resources):
       structure = get_structure_SnO2()
       launch_aiida(structure, code_string, resources)


   if __name__ == '__main__':
       code_string = 'vasp544mpi@gpu'
       resources = {'parallel_env': 'mpi*', 'tot_num_mpiprocs': 12}
       main(code_string, resources)

Once the calculation is done, with your PK, we can watch the results::

   In [1]: n = load_node(<PK>)

   In [2]: n.outputs.energies
   Out[2]: <ArrayData: uuid: 7fbf60d5-fd7c-4d45-9adb-af81e7348921 (pk: 165915)>

   In [3]: n.outputs.energies.get_arraynames()
   Out[3]: ['energy_no_entropy']

   In [4]: n.outputs.energies.get_array('energy_no_entropy')
   Out[4]: array([-39.90289213])

   In [5]: n.outputs.stress.get_arraynames()
   Out[5]: ['final']

   In [6]: n.outputs.stress.get_array('final')
   Out[6]:
   array([[ 2.65188465, -0.        ,  0.        ],
          [ 0.        ,  2.65188465,  0.        ],
          [ 0.        ,  0.        , -2.54327698]])

When we want to fully relax a crystal structure, the above script is
modified as follows:

1. Replace ``WorkflowFactory('vasp.vasp')`` by ``WorkflowFactory('vasp.relax')``
2. Remove ``IBRION`` from ``incar_dict``
3. Add the following setting::

       builder.relax = Bool(True)
       builder.force_cutoff = Float(1e-5)
       builder.steps = Int(10)
       builder.positions = Bool(True)  # Relax atomic positions
       builder.shape = Bool(True)      # Relax cell shape (alpha, beta, gamma)
       builder.volume = Bool(True)     # Relax volume
       builder.verbose = Bool(True)

The lattice parameters of the relax crystal structure is found by

::

   In [1]: n = load_node(<PK>)

   In [2]: n.outputs.structure_relaxed.cell
   Out[2]: [[4.77533981, 0.0, 0.0], [0.0, 4.77533981, 0.0], [0.0, 0.0, 3.21639965]]

There are more options for the relax workchain, e.g., running VASP
several time iteratively until convergence, which is used in the bulk
modulus example in the next section.

After the relaxation, sometimes the crystal symmetry can be slightly
broken by the VASP calculation, especially for hexagonal crystals. So
it is recommended to symmetrize the final structure if this is the case.
