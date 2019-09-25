.. _bulk_modulus:

===========================
Writing bulk modulus script
===========================

This section presents an example to calculate bulk modulus of
wurtzite-type SiC using a python script with launching several
AiiDA-VASP workchains. In the script, QueryBuilder and Group are used
to manage the workflow of this calculation.


Workflow
--------

Here one example is presented. Bulk modulus is calculated by the
following steps (the full script is attached at the end of this page):

1. Run a relax workchain to fully relax crystal structure
2. Wait until (1) finishes
3. Create two structures at fixed volumes with +/- 1% from the relaxed
   structure obtained at the step (1).
4. Submit two relax workchains to relax shapes of the crystal
   structures created at the step (3).
5. Compute bulk modulus as a post process by the formula :math:`K \simeq -V_0
   \frac{\Delta P}{\Delta  V}`

If this calculation will be done repeated and robustly, the workflow
should be written as a workchain. But as an attempt calculation as a
part of the process of daily research, writing a simple script like
following is useful by employing basic AiiDA features.


Bulk modulus calculation without using AiiDA-VASP
--------------------------------------------------

Steps 1 and 2
^^^^^^^^^^^^^

POSCAR file

::

   wurtzite-type SiC
     1.0000000000
     3.0920000000   0.0000000000   0.0000000000
    -1.5460000000   2.6777505485   0.0000000000
     0.0000000000   0.0000000000   5.0730000000
   Si    C
       2     2
   Direct
     0.3333333333   0.6666666667   0.0000000000
     0.6666666667   0.3333333333   0.5000000000
     0.3333333333   0.6666666667   0.3758220000
     0.6666666667   0.3333333333   0.8758220000

INCAR file

::

   EDIFF = 1e-08
   EDIFFG = -1e-05
   ENCUT = 500
   GGA = PS
   IALGO = 38
   IBRION = 2
   ISIF = 3
   ISMEAR = 0
   LCHARG = .FALSE.
   LREAL = .FALSE.
   LWAVE = .FALSE.
   NELM = 100
   NELMIN = 5
   NSW = 10
   PREC = Accurate
   SIGMA = 0.01

KPOINTS file

::

   # Half grid shift along c*
   0
   Gamma
               6             6             4
     0.000000000   0.000000000   0.500000000

Using this setting files, we get CONTCAR::

   SiC
      1.00000000000000
        3.0779853535726360    0.0000000000000000    0.0000000000000000
       -1.5389926767863180    2.6656135086688661    0.0000000000000000
        0.0000000000000000   -0.0000000000000000    5.0493167306164031
      Si   C
        2     2
   Direct
     0.3333333332999970  0.6666666667000030 -0.0000414569885531
     0.6666666667000030  0.3333333332999970  0.4999585430114469
     0.3333333332999970  0.6666666667000030  0.3758634569885525
     0.6666666667000030  0.3333333332999970  0.8758634569885526

     0.00000000E+00  0.00000000E+00  0.00000000E+00
     0.00000000E+00  0.00000000E+00  0.00000000E+00
     0.00000000E+00  0.00000000E+00  0.00000000E+00
     0.00000000E+00  0.00000000E+00  0.00000000E+00

Steps 3 and 4
^^^^^^^^^^^^^

We create two sets of the VASP calculation inputs. The 2nd line of
CONTCAR obtained at the step (1) is modified by the strain=0.99 (i.e.,
the 2nd line value is :math:`0.99^{1/3}` = 0.9966554934125964) and
strain=1.01 (i.e., the 2nd line value is :math:`1.01^{1/3}` =
1.0033222835420892) to change the volumes of the cells and saved as
POSCAR. INCAR is modified to have ``ISIF = 4`` to relax the cells with
keeping the volumes.

After running VASP calculations, we find the following values in the vasprun.xml's:

- strain=0.99 (volume = 41.01394436)::

       <varray name="stress" >
        <v>      22.73458454       0.00000000       0.00000000 </v>
        <v>       0.00000000      22.73458454       0.00000000 </v>
        <v>       0.00000000       0.00000000      22.73469456 </v>
       </varray>

- strain=1.01 (volume = 41.84250889)::

       <varray name="stress" >
        <v>     -21.66753480      -0.00000000      -0.00000000 </v>
        <v>       0.00000000     -21.66753480       0.00000000 </v>
        <v>       0.00000000       0.00000000     -21.66848806 </v>
       </varray>

Step 5
^^^^^^

The bulk modulus is obtained from these results as

::

   In [1]: -(41.84250889 + 41.01394436) / 2 * ((-21.66753480 * 2 - 21.66848806) / 3 - (22.73458454 * 2 + 22.73469456) / 3) / (41.84250889 - 41.01394436) / 10
   Out[1]: 222.0123695032054

So we should get bulk modulus of ~222 GPa by using AiiDA-VASP, too, as
shown below.


AiiDA-VASP script
-----------------

::

   def main(code_string, resources, group_name, sleep_seconds=60):
       group = load_group(group_name)
       structure = get_structure_SiC()
       node_relax = launch_aiida_full_relax(structure, code_string, resources,
                                            "SiC VASP calc to relax volume")
       group.add_nodes(node_relax)

       while True:
           if node_relax.is_terminated:
               break
           print("Waiting for relaxation calculation to be done.")
           sleep(sleep_seconds)

       if node_relax.is_finished_ok:
           for strain, label in zip((0.99, 1.01), ("minus", "plus")):
               structure = node_relax.outputs.structure_relaxed.clone()
               structure.set_cell(np.array(structure.cell) * strain ** (1.0 / 3))
               node = launch_aiida_relax_shape(
                   structure, code_string, resources,
                   "SiC VASP relax shape at %s volume (%f)" % (label, strain))
               group.add_nodes(node)
               print(node)
       else:
           print("Relaxation calculation failed.")


From the result of this calculation, the bulk modulus is computed by::

   import numpy as np
   from aiida.manage.configuration import load_profile
   from aiida.orm import Group, QueryBuilder
   load_profile()


   def calc_bulk_modulus(group_name):
       qb = QueryBuilder()
       qb.append(Group, filters={'label': {'==': group_name}})
       if qb.count() == 0:
           raise RuntimeError("Group %s doesn't exist." % group_name)

       stresses = []
       volumes = []
       for comment in ("minus", "plus"):
           qb = QueryBuilder()
           qb.append(Group, filters={'label': {'==': group_name}}, tag='group')
           qb.append(WorkChainNode, with_group='group',
                     filters={'label': {'ilike': '%' + comment + '%'}})
           node = qb.first()[0]
           stresses.append(np.trace(node.outputs.stress.get_array('final')) / 3)
           volumes.append(np.linalg.det(node.inputs.structure.cell))

       d_s = stresses[1] - stresses[0]
       d_v = volumes[1] - volumes[0]
       v0 = (volumes[0] + volumes[1]) / 2
       bulk_modulus = - d_s / d_v * v0

       print("Bulk modules: %f GPa" % (bulk_modulus / 10))


   if __name__ == '__main__':
       calc_bulk_modulus("Bulk modulues example")

We get the value::

   Bulk modules: 222.016084 GPa


Migration from a simple script to the WorkChain
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the calculation above, the VASP full relax workchain and the two
VASP relax workchain at constant volue are performed independently and
those nodes are just grouped. This means the workflow in ``main``
method is lost. To keep the workflow, the next challenge will be
writing a workchain of this workflow. This migration from the simple
script to the workchain will be straightforward, once we confirm the
simple script starts to work.


Full script to compute bulk modulus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   from time import sleep
   import numpy as np
   from aiida.manage.configuration import load_profile
   from aiida.orm import (
       Bool, Int, Float, Str, Code, load_group, QueryBuilder, Group,
       WorkChainNode)
   from aiida.plugins import DataFactory, WorkflowFactory
   from aiida.engine import submit
   load_profile()


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


   def launch_aiida_relax_shape(structure, code_string, resources, label):
       Dict = DataFactory('dict')
       KpointsData = DataFactory("array.kpoints")
       base_incar_dict = {
           'PREC': 'Accurate',
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

       base_config = {'code_string': code_string,
                      'potential_family': 'PBE.54',
                      'potential_mapping': {'Si': 'Si', 'C': 'C'},
                      'options': {'resources': resources,
                                  'max_wallclock_seconds': 3600 * 10}}
       base_parser_settings = {'add_energies': True,
                               'add_forces': True,
                               'add_stress': True}
       code = Code.get_from_string(base_config['code_string'])
       Workflow = WorkflowFactory('vasp.relax')
       builder = Workflow.get_builder()
       builder.code = code
       builder.parameters = Dict(dict=base_incar_dict)
       builder.structure = structure
       builder.settings = Dict(dict={'parser_settings': base_parser_settings})
       builder.potential_family = Str(base_config['potential_family'])
       builder.potential_mapping = Dict(dict=base_config['potential_mapping'])
       kpoints = KpointsData()
       kpoints.set_kpoints_mesh([6, 6, 4], offset=[0, 0, 0.5])
       builder.kpoints = kpoints
       builder.options = Dict(dict=base_config['options'])
       builder.metadata.label = label
       builder.metadata.description = label
       builder.clean_workdir = Bool(False)
       builder.relax = Bool(True)
       builder.force_cutoff = Float(1e-5)
       builder.steps = Int(10)
       builder.positions = Bool(True)
       builder.shape = Bool(True)
       builder.volume = Bool(False)
       builder.verbose = Bool(True)
       node = submit(builder)
       return node


   def launch_aiida_full_relax(structure, code_string, resources, label):
       Dict = DataFactory('dict')
       KpointsData = DataFactory("array.kpoints")
       base_incar_dict = {
           'PREC': 'Accurate',
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

       base_config = {'code_string': code_string,
                      'potential_family': 'PBE.54',
                      'potential_mapping': {'Si': 'Si', 'C': 'C'},
                      'options': {'resources': resources,
                                  'max_wallclock_seconds': 3600 * 10}}
       base_parser_settings = {'add_energies': True,
                               'add_forces': True,
                               'add_stress': True}
       code = Code.get_from_string(base_config['code_string'])
       Workflow = WorkflowFactory('vasp.relax')
       builder = Workflow.get_builder()
       builder.code = code
       builder.parameters = Dict(dict=base_incar_dict)
       builder.structure = structure
       builder.settings = Dict(dict={'parser_settings': base_parser_settings})
       builder.potential_family = Str(base_config['potential_family'])
       builder.potential_mapping = Dict(dict=base_config['potential_mapping'])
       kpoints = KpointsData()
       kpoints.set_kpoints_mesh([6, 6, 4], offset=[0, 0, 0.5])
       builder.kpoints = kpoints
       builder.options = Dict(dict=base_config['options'])
       builder.metadata.label = label
       builder.metadata.description = label
       builder.clean_workdir = Bool(False)
       builder.relax = Bool(True)
       builder.force_cutoff = Float(1e-5)
       builder.steps = Int(10)
       builder.positions = Bool(True)
       builder.shape = Bool(True)
       builder.volume = Bool(True)
       builder.convergence_on = Bool(True)
       builder.convergence_volume = Float(1e-5)
       builder.convergence_max_iterations = Int(2)
       builder.verbose = Bool(True)

       node = submit(builder)
       return node


   def main(code_string, resources, group_name, sleep_seconds=60):
       group = load_group(group_name)
       structure = get_structure_SiC()
       node_relax = launch_aiida_full_relax(structure, code_string, resources,
                                            "SiC VASP calc to relax volume")
       group.add_nodes(node_relax)

       while True:
           if node_relax.is_terminated:
               break
           print("Waiting for relaxation calculation to be done.")
           sleep(sleep_seconds)

       if node_relax.is_finished_ok:
           for strain, label in zip((0.99, 1.01), ("minus", "plus")):
               structure = node_relax.outputs.structure_relaxed.clone()
               structure.set_cell(np.array(structure.cell) * strain ** (1.0 / 3))
               node = launch_aiida_relax_shape(
                   structure, code_string, resources,
                   "SiC VASP relax shape at %s volume (%f)" % (label, strain))
               group.add_nodes(node)
               print(node)
       else:
           print("Relaxation calculation failed.")


   def calc_bulk_modulus(group_name):
       stresses = []
       volumes = []
       for label in ("minus", "plus"):
           qb = QueryBuilder()
           qb.append(Group, filters={'label': group_name}, tag='group')
           qb.append(WorkChainNode, with_group='group',
                     filters={'label': {'ilike': '%' + label + '%'}})
           node = qb.first()[0]
           stresses.append(np.trace(node.outputs.stress.get_array('final')) / 3)
           volumes.append(np.linalg.det(node.inputs.structure.cell))

       d_s = stresses[1] - stresses[0]
       d_v = volumes[1] - volumes[0]
       v0 = (volumes[0] + volumes[1]) / 2
       bulk_modulus = - d_s / d_v * v0

       print("Bulk modules: %f GPa" % (bulk_modulus / 10))


   if __name__ == '__main__':
       # code_string is chosen among the list given by 'verdi code list'
       code_string = 'vasp544mpi@gpu'

       # metadata.options.resources
       # See https://aiida.readthedocs.io/projects/aiida-core/en/latest/scheduler/index.html
       # resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 20}
       resources = {'parallel_env': 'mpi*', 'tot_num_mpiprocs': 12}

       # Here it assumes existance of the group "Bulk_modulus_SiC_test",
       # made by 'verdi group creat "Bulk_modulus_SiC_test"'.
       group_name  = "Bulk_modulus_SiC_test"
       main(code_string, resources, group_name)
       # calc_bulk_modulus(group_name)
