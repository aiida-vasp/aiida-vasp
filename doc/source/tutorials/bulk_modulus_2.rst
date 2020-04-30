exi.. _bulk_modulus_workchain:

=================================
7. Writing bulk modulus workchain
=================================

In the :ref:`previous tutorial <bulk_modulus_script>`, the initial
relaxation and the two volume restricted relaxations are performed
independently. The resulting nodes are just grouped and we execute a
calculation of the bulk modulus by fetching these nodes later.  This
means the workflow in ``main`` method is not very transparent. In fact
we would like it to follow the stepwize definition at the top as
closely as possible.

Can we do better than the scripts above? Of course. AiiDA has the
concept of a :ref:`workchains` which basically is a container for a
workflow. Not only would we like to write this workchain as modular
and reusable as possible, such that ultimately several :ref:`WorkChains` can
be cherry picked into a bigger composition of some kind of master
piece of a workchain to solve some given problem.

Let us try to preserve the workflow. The next challenge will be
writing a suitable workchain of this workflow. The migration from the
:ref:`previous script <bulk_modulus_script>` to a  will be rather
straightforward.

We will in this tutorial also show how it is possible to a plugin that
contains a workchain. This nicely demonstrate the modularity of AiiDA.

At https://github.com/atztogo/aiida-vasp-bm, the workchain
(``aiida_vasp_bm/workchains/bulkmodulus.py``) and the launch script
(``aiida_vasp_bm/example/submit_SiC.py``) shown below are obtained.

The important note about running workchain is found in `AiiDA
documentation
<https://aiida-core.readthedocs.io/en/latest/working/workflows.html#launching-work-chains>`_
or `AiiDA tutorial (see scrolling down to Warning block )
<https://aiida-tutorials.readthedocs.io/en/latest/pages/2019_SINTEF/sections/workflows.html#workchains-or-how-not-to-get-lost-if-your-computer-shuts-down-or-crashes>`_. As
written these documentations, written :ref:`workchains` have to be exposed as a
python module to be seen from AiiDA daemon. This is also
achieved by making the AiiDA plugin and installing it by ``pip``. To
make the plugin, use of AiiDA plugin cutter
(https://github.com/aiidateam/aiida-plugin-cutter) is the recommended
starting point. In this way, aiida-vasp-bm was created and this can be
installed by

::

   git clone https://github.com/atztogo/aiida-vasp-bm.git
   cd aiida-vasp-bm
   pip install -e .
   reentry scan
   verdi daemon restart


Workflow 2.0
------------

The workflow is same as the :ref:`previous
one<workflow_bulk_modulus>`, but we will now be a bit more explicit to
comply with a typical implementation of a WorkChain::

1. ``initialize`` Initialize whatever needs initializing
2. ``run_relax`` Run a relax workchain to fully relax crystal structure
3. ``create_two_structures`` Create two structures with a +/- 1% change
   in the volume from the relaxed structure obtained at step (2).
4. ``run_two_volumes`` Submit two relaxations to relax the shape of the
   structures created at step (3).
5. ``calc_bulk_modulus`` Compute bulk modulus as a post process by using the
   formula :math:`K \simeq -V_0 \frac{\Delta P}{\Delta V}`, where the
   pressure :math:`P \equiv \mathrm{Tr}(\sigma)/3` follows the VASP
   convention where :math:`\sigma` is the stress tensor.

The names, i.e. ``initialize`` above are the Python method (or
function) names in the bulk modulus WorkChain shown below. The concept
Because of WorkChain, the step (3) is executed after finishing the
step (2).

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
  workflow. This is often done by calcfuntion. For this purpose, the
  strained crystal structures are created in
  ``get_strained_structure`` and the bulk modulus is calculated in
  ``calculate_bulk_modulus`` with decorated by ``@calcfunction``, and
  these are called from the methods in WorkChain.

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
   from aiida.common.extendeddicts import AttributeDict
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
                  'account:' 'nn9995k',
                  'max_memory_kb:': 10240000,
                  'max_wallclock_seconds': 3600 * 10}

       potential_family = 'pbe'
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
       relax = AttributeDict()
       relax.perform = Bool(True)
       relax.force_cutoff = Float(1e-8)
       relax.steps = Int(10)
       relax.positions = Bool(True)
       relax.shape = Bool(True)
       relax.volume = Bool(True)
       relax.convergence_on = Bool(True)
       relax.convergence_volume = Float(1e-8)
       relax.convergence_max_iterations = Int(2)
       builder.relax = relax
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
       code_string = 'vasp@saga'
       resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 20}
       main(code_string, resources)

After running this calculation, we get the bulk modulus by

::

   In [1]: n = load_node(<PK>)

   In [2]: n.outputs.bulk_modulus.value
   Out[2]: 222.01637836634
