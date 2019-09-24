.. _bulk_modulus:

=================================
Example: Bulk modulus calculation
=================================

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


When we want to relax a crystal structure, the above script is
modified as follows:

1. ``WorkflowFactory('vasp.relax')``
2. Remove ``IBRION`` from ``incar_dict``
3. Add the following setting::

       builder.relax = Bool(True)
       builder.force_cutoff = Float(1e-5)
       builder.convergence_on = Bool(True)
       builder.convergence_volume = Float(1e-5)
       builder.convergence_max_iterations = Int(10)
       builder.relax_parameters = Dict(dict={'EDIFFG': -1e-5,
                                             'IBRION': 2,
                                             'NSW': 10,
                                             'ISIF': 3})
       builder.verbose = Bool(True)

After the relaxation, somethimes the crystal symmetry can be slightly
broken by the VASP calculation, especially for hexagonal crystals. It
is recommended to symmetrize the final structure if this is minded.


Use of Group and QueryBuilder of AiiDA
---------------------------------------

Once we start daily use of AiiDA to run VASP calculations, we will
meet the problem how to remember the location of results. We are
familier with handling files in directories/folders on conventional
file system, but the data in AiiDA are stored in the database.

Group
^^^^^

The initial easiest choice to get similar feeling to directories is
the use of Group. We make groups (``verdi group create``) and put
workchain nodes into them. The details are found at the `official
documentation
<https://aiida-core.readthedocs.io/en/latest/working_with_aiida/groups.html>`_.

QueryBuilder
^^^^^^^^^^^^

The next step is the use of QueryBuilder. This offers to search nodes
with given hints such as label, created time, and node type. This
definitely provides flexible search of data. For the begginers, it may
be painful to use it, however we have to learn how to use it for our
vigorous life. The official documentation for QueryBuilder is found
`here
<https://aiida-core.readthedocs.io/en/latest/working_with_aiida/index.html#querying-data>`_,
but the `tutorial material <https://aiida-tutorials.readthedocs.io/en/tutorial_sintef/pages/2019_SINTEF/sections/querybuilder.html>`_ would give a better catch.

Running VASP calculations after VASP relax calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here one example is presented. Bulk modulus is calculated by the
following steps (the full script is attached at the end of this page):

1. Run VASP relax workchain
2. Wait until (1) finishes
3. Submit two VASP calculations with ``ISIF=4`` to get stresses at volume
   ratios of 0.99 and 1.01 with respect to the relaxed structure.
4. Compute bulk modulus as a post process by the formula :math:`K \simeq -V_0
   \frac{\Delta P}{\Delta  V}`

If this calculation will be done repeated and robustly, the workflow
should be written as a workchain. But as an attempt calculation as a
part of the process of daily research, writing a simple script like
following is useful by employing basic AiiDA features.

::

   def main(code_string, resources,
            group_name="Bulk modulues example",
            sleep_seconds=60):
       qb = QueryBuilder()
       qb.append(Group, filters={'label': {'==': group_name}})
       if qb.count() == 0:
           group = Group(label=group_name)
           group.store()
           print("Group %s was created." % group_name)
       else:
           group = load_group(group_name)
       structure = get_structure_AlN()
       node_relax = launch_aiida_relax(structure, code_string, resources,
                                       label="AlN VASP calc to relax volume")
       group.add_nodes(node_relax)

       while True:
           if node_relax.is_terminated:
               break
           print("Waiting for relaxation calculation to be done.")
           sleep(sleep_seconds)

       if node_relax.is_finished_ok:
           for strain, comment in zip((0.99, 1.01), ('minus', 'plus')):
               structure = node_relax.outputs.structure_relaxed.clone()
               structure.set_cell(np.array(structure.cell) * strain ** (1.0 / 3))
               label = "AlN VASP calc at %s volume (%f)" % (comment, strain)
               node = launch_aiida(structure, code_string, resources, label=label)
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

   Bulk modules: 201.982655 GPa


From a simple script to workchain
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the calculation above, the VASP relax calculation and the two VASP
calculations are independently calculated and just grouped. This means
the workflow is lost. The next challenge will be writing the workchain
of this workflow.

Migration will be straightforward, once this simple script starts to
work and how to design and write workchains are understood.


Full script to compute bulk modulus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   from time import sleep
   import numpy as np
   from aiida.manage.configuration import load_profile
   from aiida.orm import (
       Bool, Str, Code, Int, Float, load_group, Group,
       QueryBuilder, WorkChainNode)
   from aiida.plugins import DataFactory, WorkflowFactory
   from aiida.engine import submit
   load_profile()


   def launch_aiida(structure, code_string, resources,
                    label="AlN VASP calculation"):
       Dict = DataFactory('dict')
       KpointsData = DataFactory("array.kpoints")

       incar_dict = {'PREC': 'Accurate',
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
                     'IBRION': 2,
                     'NSW': 10,
                     'ISIF': 4,
                     'EDIFFG': -1e-8}

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


   def launch_aiida_relax(structure, code_string, resources,
                          label="AlN VASP relax calculation"):
       Dict = DataFactory('dict')
       KpointsData = DataFactory("array.kpoints")

       incar_dict = {
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
       Workflow = WorkflowFactory('vasp.relax')
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
       builder.relax = Bool(True)
       builder.force_cutoff = Float(1e-5)
       builder.convergence_on = Bool(True)
       builder.convergence_volume = Float(1e-5)
       builder.convergence_max_iterations = Int(10)
       builder.relax_parameters = Dict(dict={'IBRION': 2,
                                             'NSW': 10,
                                             'ISIF': 3,
                                             'EDIFFG': -1e-8})
       builder.verbose = Bool(True)

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


   def main(code_string, resources,
            group_name="Bulk modulues example",
            sleep_seconds=60):
       qb = QueryBuilder()
       qb.append(Group, filters={'label': {'==': group_name}})
       if qb.count() == 0:
           group = Group(label=group_name)
           group.store()
           print("Group %s was created." % group_name)
       else:
           group = load_group(group_name)
       structure = get_structure_AlN()
       node_relax = launch_aiida_relax(structure, code_string, resources,
                                       label="AlN VASP calc to relax volume")
       group.add_nodes(node_relax)

       while True:
           if node_relax.is_terminated:
               break
           print("Waiting for relaxation calculation to be done.")
           sleep(sleep_seconds)

       if node_relax.is_finished_ok:
           for strain, comment in zip((0.99, 1.01), ('minus', 'plus')):
               structure = node_relax.outputs.structure_relaxed.clone()
               structure.set_cell(np.array(structure.cell) * strain ** (1.0 / 3))
               label = "AlN VASP calc at %s volume (%f)" % (comment, strain)
               node = launch_aiida(structure, code_string, resources, label=label)
               group.add_nodes(node)
               print(node)
       else:
           print("Relaxation calculation failed.")


   if __name__ == '__main__':
       code_string = 'vasp544mpi@gpu'
       resources = {'parallel_env': 'mpi*', 'tot_num_mpiprocs': 12}
       main(code_string, resources)
