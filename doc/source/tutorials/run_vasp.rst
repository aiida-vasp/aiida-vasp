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

The usage of spglib is found at https://atztogo.github.io/spglib/.

For example, POSCAR of wurtzite-AlN structure can be obtained from
the Materials project database,
https://materialsproject.org/materials/mp-661/#, which is::

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

This is already symmetrized but we can symmetrize more if we
want more number of digits for :math:`\sqrt{3}` and
:math:`\frac{1}{3}` that appear in hexagonal crystals.
An ad-hoc script to symmetrize this may be::

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

A typical script to launch a VASP caluculation is something like::

   import numpy as np
   from aiida.manage.configuration import load_profile
   from aiida.orm import Bool, Str, Code
   from aiida.plugins import DataFactory, WorkflowFactory
   from aiida.engine import submit
   load_profile()


   def launch_aiida(structure, code_string, resources,
                    label="AlN VASP calculation"):
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
                  'max_wallclock_seconds': 3600 * 10}

       potential_family = 'PBE.54'
       potential_mapping = {'Al': 'Al', 'N': 'N'}

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


   def get_structure_AlN():
       """Set up AlN primitive cell

        Al N
          1.0
            3.1100000000000000    0.0000000000000000    0.0000000000000000
           -1.5550000000000000    2.6933390057696038    0.0000000000000000
            0.0000000000000000    0.0000000000000000    4.9800000000000000
        Al N
          2   2
       Direct
          0.3333333333333333  0.6666666666666665  0.0000000000000000
          0.6666666666666667  0.3333333333333333  0.5000000000000000
          0.3333333333333333  0.6666666666666665  0.6190000000000000
          0.6666666666666667  0.3333333333333333  0.1190000000000000

       """

       StructureData = DataFactory('structure')
       a = 3.11
       c = 4.98
       lattice = [[a, 0, 0],
                  [-a / 2, a / 2 * np.sqrt(3), 0],
                  [0, 0, c]]
       structure = StructureData(cell=lattice)
       for pos_direct, symbol in zip(
               ([1. / 3, 2. / 3, 0],
                [2. / 3, 1. / 3, 0.5],
                [1. / 3, 2. / 3, 0.619],
                [2. / 3, 1. / 3, 0.119]), ('Al', 'Al', 'N', 'N')):
           pos_cartesian = np.dot(pos_direct, lattice)
           structure.append_atom(position=pos_cartesian, symbols=symbol)
       return structure


   def main(code_string, resources):
       structure = get_structure_AlN()
       launch_aiida(structure, code_string, resources)


   if __name__ == '__main__':
       code_string = 'vasp544mpi@gpu'
       resources = {'parallel_env': 'mpi*', 'tot_num_mpiprocs': 12}
       main(code_string, resources)


When we want to fully relax a crystal structure, the above script is
modified as follows:

1. ``WorkflowFactory('vasp.relax')``
2. Remove ``IBRION`` from ``incar_dict``
3. Add the following setting::

       builder.relax = Bool(True)
       builder.force_cutoff = Float(1e-5)
       builder.steps = Int(10)
       builder.positions = Bool(True)  # Relax atomic positions
       builder.shape = Bool(True)      # Relax cell shape (alpha, beta, gamma)
       builder.volume = Bool(True)     # Relax volume
       builder.verbose = Bool(True)

After the relaxation, sometimes the crystal symmetry can be slightly
broken by the VASP calculation, especially for hexagonal crystals. It
is recommended to symmetrize the final structure if this is minded.
