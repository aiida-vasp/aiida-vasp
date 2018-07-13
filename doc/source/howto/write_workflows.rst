How to write an AiiDA-VASP workflow
===================================

Before we start
---------------

This document describes how to creae workflows for VASP using AiiDA and the AiiDA-VASP plugin. It assumes the reader to be familiar with the `AiiDA docs`_ on how to write workflows as well as with `VASP`_

.. _AiiDA docs: https://aiida-core.readthedocs.io/en/stable/work/index.html
.. _VASP: https://www.vasp.at/index.php/documentation


Modularity
----------

AiiDA workflows are modular and can be nested indefinitely. This means it makes sense to build a hierarchy of reusable workflows, where the more complex workflows reuse a set of simple ones. This example will show how to build a fairly simple one, but the process of building more complex ones is much the same.

AiiDA-VASP already provides a VaspBaseWf (accessible via ``WorkflowFactory('vasp.base')``), which should be used to run VASP from inside other workflows, it already provides some input validation and will be extended to provide automatic restarts and error handling in the future. It also supports both AiiDA versions before and after v1.0.0 by calling the ``VaspCalculation`` in the appropriate way in both cases.

Reusable skeleton code
----------------------

From the AiiDA documentation you should be familiar with this code. It runs a VASP calculation with most parameters chosen by default and some defaults overridable::

   from aiida.work.workchain import WorkChain, submit
   from aiida.work.db_types import Str, Int
   from aiida.orm import WorkchainFactory
   from aiida.common.extendeddicts import AttributeDict

   from aiida_vasp.utils.aiida_utils import get_data_node, get_data_class

   class ExampleVaspWorkflow(WorkChain):

      @classmethod
      def define(cls, spec):
         super(ExampleVaspWorkflow, cls).define(spec)
         spec.input('code', valid_type=Code)
         spec.input('structure', valid_type=(get_data_class('structure'), get_data_class('cif')))
         spec.input('potcar_family', valid_type=Str)
         spec.input('potcar_mapping', valid_type=get_data_class('parameter'))
         spec.input('incar_overrides', valid_type=get_data_class('parameter'), required=False)
         spec.input('machines', valid_type=Int, required=False)
         spec.input('mpi_procs', valid_type=Int, required=False)
         spec.input('queue_name', valid_type=Str, required=False)

         cls.outline(
            cls.generate_missing_inputs,
            cls.run_calculation,
            cls.results
         )

         cls.output(output_parameters)

      def generate_missing_inputs(self):
         """Here we create the inputs required to rune the 'vasp.base' workflow and store them in the context."""
         self.ctx.inputs = AttributeDict()
         self.ctx.inputs.code = self.inputs.code
         self.ctx.inputs.structure = self.inputs.structure
         self.ctx.inputs.potcar_family = self.inputs.potcar_family
         self.ctx.inputs.potcar_mapping = self.inputs.potcar_mapping
         self.ctx.inputs.kpoints = self._generate_kpoints()
         self.ctx.inputs.incar = self._generate_incar(self.inputs.incar_overrides.get_dict())
         self.ctx.inputs.incar = self._generate_options(
            self.inputs.get('machines', None),
            self.inputs.get('mpi_procs', None),
            self.inputs.get('queue_name', None)
         )
         ## Follow the same pattern to use additional inputs like wave functions or charge densities

      def run_calculation(self):
         running = submit(WorkflowFactory('vasp.base'), **self.ctx.inputs)
         return ToContext(vasp_run=running)

      def results(self):
         if 'output_structure' not in self.ctx.vasp_run.out:
                self.report('"structure" was not an output of {} (pk: {})'.format(self.ctx.vasp_run, self.ctx.vasp_run.pk))
         self.out('structure', self.ctx.vasp_run.out.output_structure)

      def _generate_options(self, machines=Int(1), mpi_procs=Int(1), queue_name=Str('batch')):
         options = AttributeDict()
         options.resources = {'num_machines': machines.value, 'num_mpiprocs_per_machine': mpi_procs.value}
         options.queue_name = queue_name.value
         return get_data_node('parameter', dict=options)
      
      def _generate_kpoints(self):
         """How to generate KpointData nodes is documented in AiiDA's documentation."""
         kpoints = get_data_class('array.kpoints')
         ## ... generate and set a kpoints path or mesh
         return kpoints

      def _generate_incar(overrides):
         """Define the INCAR parameters, keys can be lower case"""
         incar_dict = {
            system: "Example Workflow system",
            isif: 4,
            ## Set defaults here
         }
         incar_dict.update(overrides)  ## Apply user provided optional overrides
         incar = get_data_node('parameter', dict=incar_dict)
         return incar

This example uses the ``vasp.base`` workflow to run a single VASP calculation with defaults. Higher complexity can be achieved using WorkChain control flow features like conditionals, loops, etc, described in the AiiDA documentation linked above.

Determine the inputs and outputs
--------------------------------

One of the first questions in designing a workflow should be which inputs will be required and what outputs should be generated. For example: a relaxation workflow will obviously require a code and a structure, as well as some way of describing which POTCAR potentials to use.

 * code: tells us which VASP executable we will run on what machine
 * structure: describes the structure to be relaxed

It might provide defaults for everything else and provide optional inputs for changing them, or the choice can be made to require some of the other parameters (for example k-point mesh density, compute resources, etc)

Determine the required steps
----------------------------

It is helpful to sketch out a flow diagram before approaching writing a workflow. How to translate such a flow diagram into a ``WorkChain`` outline should be obvious from AiiDA's documentation (linked above).

As a (simple) example: Calculating a band structure requires an (optional) relaxation of the structure, followed by an SCF run, and a band structure run using the charge densities from the SCF run.

::

   input structure -> relax -> output structure -> scf run -> chgcar -> band structure run -> band structure
                                       |                                 ^
                                       +---------------------------------+


.. _howto/base_wf/reference
Detailed usage of VaspBaseWf
----------------------------

A note about compatibility: WorkChains provide a handy pattern for interactively building input sets both under AiiDA < 1.0.0 as from AiiDA 1.0.0a1 onwards. They are very similar but different enough to recommend using a python dictionary or ``aiida.common.extendeddicts.AttributeDict`` instead in scripts where compatibility for both should be achieved.

Required inputs
^^^^^^^^^^^^^^^

The VaspBaseWf requires a number of inputs, these comprise the minimum set of information to run a VASP calculation from AiiDA.

 * ``code``: an AiiDA ``Code``, describes the VASP executable and holds a reference to the ``Computer`` instance on which it lives.
 * ``structure``: an AiiDA ``StructureData`` or ``CifData``, describes the structure on which VASP is to be run.
 * ``kpoints``: an AiiDA ``KpointsData`` instance, describing the kpoints mesh or path.
 * ``potcar_family``: an AiiDA ``Str``, the name given to a set of uploaded POTCAR files.
 * ``potcar_mapping``: an AiiDA ``ParameterData``, containing an entry for at least every kind name in the ``structure`` input with the full name of the POTCAR from the ``potcar_family``. Example: ``{'In1': 'In_d', 'In2': 'In_h'}``.
 * ``incar``: an AiiDA ``ParameterData`` instance, containing key/value pairs that get written to INCAR as ``KEY = VALUE``, keys can be lower case and builtin python types should be used for values.
 * ``options``, an AiiDA ``ParameterData`` instance, containing at least the keys ``resources`` and ``queue_name``. More information about calculation options is available in the AiiDA documentation.

Optional inputs
^^^^^^^^^^^^^^^

Optional inputs are not required and can be used to change aspects of the VASP run:

 * ``wavecar``: an instance of ``aiida_vasp.data.wavefun.WavefunData`` (factory string: ``vasp.wavefun``). Used to pass Wavefunctions from a previous run to a follow up calculation.
 * ``chgcar``: an instance of ``aiida_vasp.data.chargedensity.ChargedensityData`` (factory string: ``vasp.chargedensity``. Used to pass charge densities calculated in a previous run.
 * ``settings``: ``ParameterData``, contains additional settings for AiiDA-side aspects of the VASP calculation, like additional files to retrieve, optional quantities to be parsed, etc.

Outputs
^^^^^^^

The outputs, if no additional ones are requested using the ``settings`` input, are:

 * ``output_parameters``: ``ParameterData``, scalar and low dimensional vector quantities, like energies, forces, etc, parsed from OUTCAR and vasprun.xml
 * ``output_structure``: ``StructureData``, what VASP outputs in CONTCAR
 * ``retrieved``: ``FolderData`` containing the retrieved files
 * ``remote_folder``: ``RemoteFolderData`` containing information about the remote work folder in which VASP was run
