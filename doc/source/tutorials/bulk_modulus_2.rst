.. _bulk_modulus_workchain:

=================================
5. Writing bulk modulus workchain
=================================

In the previous :ref:`bulk_modulus_script`, the initial relaxation and the two volume restricted
relaxations are performed independently. The resulting nodes are just grouped and
we execute a calculation of the bulk modulus by fetching these nodes later.
This means the workflow in ``main``
method is not very transparent. In fact we would like it to follow the stepwize
definition at the top as closely as possible.

Can we do better than the scripts above? Of course. AiiDA has the concept of a WorkChain which
basically is a container for a workflow. Not only would we like to write this WorkChain
as modular and reusable as possible, such that
ultimately several WorkChains can be cherry picked into a bigger composition of some kind
of master piece of a workchain to solve some given problem.

Let us try to preserve the workflow. The next challenge will be
writing a suitable WorkChain of this workflow. This migration from the
script to the workchain will be straightforward. But first, please study and try the script
below, which is the complete script to calculate the bulk modulus.

This section again presents an example to calculate bulk modulus of
wurtzite-type SiC by writing a workchain. We migrate what we did in the
:ref:`last section <bulk_modulus>` to a WorkChain.

At https://github.com/atztogo/aiida-vasp-bm, the WorkChain
(``aiida_vasp_bm/workchains/bulkmodulus.py``) and the launch script
(``aiida_vasp_bm/example/submit_SiC.py``) shown below are obtained.


Workflow
--------

The workflow is same as the previous one, but the steps are put in
WorkChain:

1. [``run_relax``] Run a relax workchain to fully relax crystal structure
2. [``create_two_structures``] Create two structures at fixed volumes
   with +/- 1% from the relaxed structure obtained at the step (1).
3. [``run_two_volumes``] Submit two relax workchains to relax shapes of the crystal
   structures created at the step (2).
4. [``run_two_volumes``] Compute bulk modulus as a post process by the
   formula :math:`K \simeq -V_0 \frac{\Delta P}{\Delta V}`, where the
   pressure :math:`P \equiv \mathrm{Tr}(\sigma)/3` following the VASP
   convention with :math:`\sigma` as the stress tensor.

The names just after the numbers in the above list are the method
names in the bulk modulus workchain shown below. Because of WorkChain,
the step (2) is executed after finishing the step (1).

Implementation of Workchain
---------------------------

- This workchain is written re-using specification of inputs of the
  relax workchain in the AiiDA-VASP plugin
  (``spec.expose_inputs(RelaxWorkChain)``).
- One output is defined for bulk modulus in GPa as a float value
  (``spec.output('bulk_modulus', valid_type=Float)``.)
- The workflow is described as::

    spec.outline(
        cls.initialize,
        cls.run_relax,
        cls.create_two_structures,
        cls.run_two_volumes,
        cls.calc_bulk_modulus,
    )

- Data created inside WorkChain should be connected to the
  workflow. This is often done by calcfuntion. For this purpose, the strained crystal
  structures are created in ``get_strained_structure`` and the bulk
  modulus is calculated in ``calculate_bulk_modulus`` with decorated
  by ``@calcfunction``, and these are called from the methods in WorkChain.

::

   import numpy as np
   from aiida.orm import Bool
   from aiida.plugins import DataFactory, WorkflowFactory
   from aiida.engine import WorkChain, calcfunction
   from aiida.common.extendeddicts import AttributeDict

   Float = DataFactory('float')


   @calcfunction
   def get_strained_structure(structure, strain):
       new_structure = structure.clone()
       new_structure.set_cell(
           np.array(new_structure.cell) * strain.value ** (1.0 / 3))
       return new_structure


   @calcfunction
   def calculate_bulk_modulus(stress_minus, stress_plus,
                              structure_minus, structure_plus):
       stresses = []
       volumes = []
       for stress in (stress_minus, stress_plus):
           stresses.append(np.trace(stress.get_array('final')) / 3)
       for structure in (structure_minus, structure_plus):
           volume = np.linalg.det(structure.cell)
           volumes.append(volume)
       d_s = stresses[1] - stresses[0]
       d_v = volumes[1] - volumes[0]
       v0 = (volumes[0] + volumes[1]) / 2
       bulk_modulus = - d_s / d_v * v0 / 10  # GPa
       return Float(bulk_modulus)


   class BulkModulusWorkChain(WorkChain):
       """WorkChain to compute bulk modulus using VASP."""

       _next_workchain_string = 'vasp.relax'
       _next_workchain = WorkflowFactory(_next_workchain_string)

       @classmethod
       def define(cls, spec):
           super(BulkModulusWorkChain, cls).define(spec)
           spec.expose_inputs(cls._next_workchain)
           spec.outline(
               cls.initialize,
               cls.run_relax,
               cls.create_two_structures,
               cls.run_two_volumes,
               cls.calc_bulk_modulus,
           )
           spec.output('bulk_modulus', valid_type=Float)

       def initialize(self):
           self.report("initialize")
           self.ctx.inputs = AttributeDict()
           self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))

       def run_relax(self):
           self.report("run_relax")
           Workflow = WorkflowFactory('vasp.relax')
           builder = Workflow.get_builder()
           for key in self.ctx.inputs:
               builder[key] = self.ctx.inputs[key]
           if 'label' in self.ctx.inputs.metadata:
               label = self.ctx.inputs.metadata['label'] + " relax"
               builder.metadata['label'] = label
           if 'description' in self.ctx.inputs.metadata:
               description = self.ctx.inputs.metadata['description'] + " relax"
               builder.metadata['description'] = description
           future = self.submit(builder)
           self.to_context(**{'relax': future})

       def create_two_structures(self):
           self.report("create_two_structures")
           for strain, name in zip((0.99, 1.01), ('minus', 'plus')):
               structure = get_strained_structure(
                   self.ctx['relax'].outputs.structure_relaxed, Float(strain))
               structure.label = name
               self.ctx['structure_%s' % name] = structure

       def run_two_volumes(self):
           self.report("run_two_volumes")
           for strain, future_name in zip((0.99, 1.01), ('minus', 'plus')):
               Workflow = WorkflowFactory('vasp.relax')
               builder = Workflow.get_builder()
               for key in self.ctx.inputs:
                   builder[key] = self.ctx.inputs[key]
               if 'label' in self.ctx.inputs.metadata:
                   label = self.ctx.inputs.metadata['label'] + " " + future_name
                   builder.metadata['label'] = label
               if 'description' in self.ctx.inputs.metadata:
                   description = self.ctx.inputs.metadata['description']
                   description += " " + future_name
                   builder.metadata['description'] = description
               builder.structure = self.ctx['structure_%s' % future_name]
               builder.force_cutoff = Float(1e-8)
               builder.positions = Bool(True)
               builder.shape = Bool(True)
               builder.volume = Bool(False)
               builder.convergence_on = Bool(False)
               future = self.submit(builder)
               self.to_context(**{future_name: future})

       def calc_bulk_modulus(self):
           self.report("calc_bulk_modulus")
           bulk_modulus = calculate_bulk_modulus(
               self.ctx['minus'].outputs.stress,
               self.ctx['plus'].outputs.stress,
               self.ctx['minus'].inputs.structure,
               self.ctx['plus'].inputs.structure)
           bulk_modulus.label = "Bulk modulus in GPa"
           self.out('bulk_modulus', bulk_modulus)
           self.report('finish bulk modulus calculation')


Launch script
-------------

::

   import numpy as np
   from aiida.manage.configuration import load_profile
   from aiida.orm import Bool, Str, Code, Int, Float, WorkChainNode, QueryBuilder, Group
   from aiida.plugins import DataFactory, WorkflowFactory
   from aiida.engine import submit

   load_profile()

   Dict = DataFactory('dict')
   KpointsData = DataFactory("array.kpoints")


   def launch_aiida_bulk_modulus(structure, code_string, resources,
                                 label="VASP bulk modulus calculation"):
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
       potential_mapping = {'Si': 'Si', 'C': 'C'}

       parser_settings = {'add_energies': True,
                          'add_forces': True,
                          'add_stress': True}

       code = Code.get_from_string(code_string)
       Workflow = WorkflowFactory('vasp_bm.bulkmodulus')
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
       builder.force_cutoff = Float(1e-8)
       builder.steps = Int(10)
       builder.positions = Bool(True)
       builder.shape = Bool(True)
       builder.volume = Bool(True)
       builder.convergence_on = Bool(True)
       builder.convergence_volume = Float(1e-8)
       builder.convergence_max_iterations = Int(2)
       builder.verbose = Bool(True)

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
       node = launch_aiida_bulk_modulus(structure, code_string, resources,
                                        label="SiC VASP bulk modulus calculation")
       print(node)


   if __name__ == '__main__':
       code_string = 'vasp544mpi@gpu'
       resources = {'parallel_env': 'mpi*', 'tot_num_mpiprocs': 12}
       main(code_string, resources)

After running this calculation, we get the bulk modulus by

::

   In [1]: n = load_node(<PK>)

   In [2]: n.outputs.bulk_modulus.value
   Out[2]: 222.01637836634
