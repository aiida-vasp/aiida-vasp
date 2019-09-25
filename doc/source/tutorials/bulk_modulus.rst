.. _bulk_modulus:

===========================
Writing bulk modulus script
===========================

This section presents an example to calculate bulk modulus of
rutile-type SnO2 using a python script with launching several
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

::

   def main(code_string, resources, group_name, sleep_seconds=60):
       group = load_group(group_name)
       structure = get_structure_SnO2()
       node_relax = launch_aiida_full_relax(structure, code_string, resources,
                                            "SnO2 VASP calc to relax volume")
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
                   "SnO2 VASP relax shape at %s volume (%f)" % (label, strain))
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

   Bulk modules: 193.577285 GPa


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
                      'kpoints_density': 0.5,  # k-point density,
                      'potential_family': 'PBE.54',
                      'potential_mapping': {'Sn': 'Sn', 'O': 'O'},
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
       kpoints.set_kpoints_mesh([4, 4, 6], offset=[0.5, 0.5, 0.5])
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
                      'kpoints_density': 0.5,  # k-point density,
                      'potential_family': 'PBE.54',
                      'potential_mapping': {'Sn': 'Sn', 'O': 'O'},
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
       kpoints.set_kpoints_mesh([4, 4, 6], offset=[0.5, 0.5, 0.5])
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


   def main(code_string, resources, group_name, sleep_seconds=60):
       group = load_group(group_name)
       structure = get_structure_SnO2()
       node_relax = launch_aiida_full_relax(structure, code_string, resources,
                                            "SnO2 VASP calc to relax volume")
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
                   "SnO2 VASP relax shape at %s volume (%f)" % (label, strain))
               group.add_nodes(node)
               print(node)
       else:
           print("Relaxation calculation failed.")


   if __name__ == '__main__':
       # code_string is chosen among the list given by 'verdi code list'
       code_string = 'vasp544mpi@gpu'

       # metadata.options.resources
       # See https://aiida.readthedocs.io/projects/aiida-core/en/latest/scheduler/index.html
       # resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 20}
       resources = {'parallel_env': 'mpi*', 'tot_num_mpiprocs': 12}

       # Here it assumes existance of the group "Bulk_modulus_SnO2_test",
       # made by 'verdi group creat "Bulk_modulus_SnO2_test"'.
       group_name  = "Bulk_modulus_SnO_test"
       # main(code_string, resources, group_name)
       calc_bulk_modulus(group_name)
