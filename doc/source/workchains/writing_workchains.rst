.. _writing_workflows:

Writing workchains
==================

Before we start
---------------

This document describes briefly how to write :ref:`workchains` that can be used as a standalone :ref:`workflows` or composed into a more extensive one. It assumes the reader to be familiar with the `AiiDA documentation`_ on how to construct a `Workchain`_ as well as how `VASP`_ operates.

Modularity
----------

:ref:`workchains` are modular and can be nested indefinitely. This means it makes sense to build a hierarchy of reusable :ref:`workchains`, where the more complex :ref:`workflows` reuse a set of fundamental ones. This example will show how to build a basic workchain, but the process of building more complex ones is similar.

Reusable skeleton code
----------------------

From the `AiiDA documentation`_ you should be familiar with this code. It runs the :ref:`VASP workchain` with most input parameters inherited from the definition specified in :ref:`VASP workchain`, but we override the ``structure``. The reason we want to do this is to demonstrate the simple example of building inputs on the fly in workchains such that it is possible for a user to call the workchain we here construct with minimal input parameters. Say that we always want to do calculations on the same structure (this could of course be a list of structures or the whole set of structures in say the Materials Project database). Since :ref:`VASP workchain` requires a ``structure`` as an input, we will in the example below override that such that we in the example workchain specify the structure. Consequently, then the user wants to call the example workchain, they do not have to supply a structure, a default is chosen if they do not::

   from aiida.engine import WorkChain
   from aiida.orm import Str, Int
   from aiida.plugins import WorkflowFactory, DataFactory
   from aiida.common.extendeddicts import AttributeDict

   class ExampleVaspWorkchain(WorkChain):

      _next_workchain = WorkflowFactory('vasp.vasp')
      
      @classmethod
      def define(cls, spec):
         super(ExampleVaspWorkchain, cls).define(spec)
	 spec.expose_inputs(cls._next_workchain, exclude=['structure'])
	 spec.input('structure', valid_type=(get_data_class('structure'), get_data_class('cif')))
	 cls.outline(
            cls.generate_inputs,
	    cls.generate_structure,
            cls.run_next_workchain,
            cls.verify_next_workchain,
	    cls.results
         )
	 spec.expose_outputs(cls._next_workchain)
	 spec.exit_code(420, 'ERROR_NO_CALLED_WORKCHAIN', message='no called workchain detected')

      def generate_missing_inputs(self):
         """Here we create the inputs required to run the 'vasp.vasp' workchain and store them in the context."""
         self.ctx.inputs = AttributeDict()
         self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))

      def generate_structure(self):
         """Here we generate the structure if it is missing from the input."""
	 try:
	     self.ctx.inputs.structure = self.inputs.structure
	 except AttributeError:
	     # Generate an example silicon structure
	     structure_class = DataFactory('structure')
	     alat = 5.4
	     structure = structure_class(cell=numpy.array([[.5, 0, .5], [.5, .5, 0], [0, .5, .5]]) * alat)
	     structure.append_atom(position=numpy.array([0.0, 0.0, 0.0]) * alat, symbols='Si')
	     structure.append_atom(position=numpy.array([.25, .25, .25]) * alat, symbols='Si')
	     self.ctx.inputs.structure = structure
	 
      def run_next_workchain(self):
         running = self.submit(self._next_workchain, **self.ctx.inputs)
         return self.to_context(workchains=running)

      def verify_next_workchain(self):
         """Make sure we attach all results coming from next_workchain to this workchain."""
	 try:
	     workchain = self.ctx.workchains[-1]
	 except IndexError:
	     self.report("Could not find the next_workchain.")
	     return self.exit_codes.ERROR_NO_CALLED_WORKCHAIN

      def results(self):
          """Attach all outputs from next_workchain to this workchain."""
	  workchain = self.ctx.workchains[-1]
	  self.out_many(self.exposed_ouputs(workchain, self._next_workchain))

This example uses the :ref:`VASP workchain` to run a single `VASP`_ calculation with its defaults. Please also consult the example files in the ``examples`` folder, which calls the bundled workchains.

Determine the inputs and outputs
--------------------------------

One of the first questions in designing a workchain should be which inputs will be required and what outputs should be generated. A workchain might provide defaults for everything and work as a passthrough, it might set up all inputs or outputs, or only parts of them. As a user writing new workchains, one should thus first be concerned about defining these and writing the ``spec.input``, ``spec.output``.


Determine the required steps
----------------------------

It is helpful to sketch out a flow diagram before approaching writing a workchain. How to translate such a flow diagram into a ``cls.outline`` should be obvious from the `AiiDA documentation`_. One should take care on trying to factor out components and avoiding to write very large workchains realize a workflow. By segmenting the problem, its steps, inputs and outputs one ensures a greater opportunity to reuse the workchain in other workflows.


Detailed usage of VaspWorkChain
-------------------------------

A note about compatibility: WorkChains provide a handy pattern for interactively building input sets both under AiiDA < 1.0.0 as from AiiDA 1.0.0a1 onwards. They are very similar but different enough to recommend using a python dictionary or ``aiida.common.extendeddicts.AttributeDict`` instead in scripts where compatibility for both should be achieved.

Required inputs
^^^^^^^^^^^^^^^

The VaspWorkChain requires a number of inputs, these comprise the minimum set of information to run a VASP calculation from AiiDA.

 * ``code``: an AiiDA ``Code``, describes the VASP executable and holds a reference to the ``Computer`` instance on which it lives.
 * ``structure``: an AiiDA ``StructureData`` or ``CifData``, describes the structure on which VASP is to be run.
 * ``kpoints``: an AiiDA ``KpointsData`` instance, describing the kpoints mesh or path.
 * ``potential_family``: an AiiDA ``Str``, the name given to a set of uploaded POTCAR files.
 * ``potential_mapping``: an AiiDA ``Dict``, containing an entry for at least every kind name in the ``structure`` input with the full name of the POTCAR from the ``potential_family``. Example: ``{'In1': 'In_d', 'In2': 'In_h'}``.
 * ``incar``: an AiiDA ``Dict`` instance, containing key/value pairs that get written to INCAR as ``KEY = VALUE``, keys can be lower case and builtin python types should be used for values.
 * ``options``, an AiiDA ``Dict`` instance, containing at least the keys ``resources`` and ``queue_name``. More information about calculation options is available in the AiiDA documentation.

Optional inputs
^^^^^^^^^^^^^^^

Optional inputs are not required and can be used to change aspects of the VASP run:

 * ``wavecar``: an instance of ``aiida_vasp.data.wavefun.WavefunData`` (factory string: ``vasp.wavefun``). Used to pass Wavefunctions from a previous run to a follow up calculation.
 * ``chgcar``: an instance of ``aiida_vasp.data.chargedensity.ChargedensityData`` (factory string: ``vasp.chargedensity``. Used to pass charge densities calculated in a previous run.
 * ``settings``: ``Dict``, contains additional settings for AiiDA-side aspects of the VASP calculation, like additional files to retrieve, optional quantities to be parsed, etc.

Outputs
^^^^^^^

The outputs, if no additional ones are requested using the ``settings`` input, are:

 * ``parameters``: ``Dict``, scalar and low dimensional vector quantities, like energies, forces, etc, parsed from OUTCAR and vasprun.xml
 * ``structure``: ``StructureData``, what VASP outputs in CONTCAR
 * ``retrieved``: ``FolderData`` containing the retrieved files
 * ``remote_folder``: ``RemoteFolderData`` containing information about the remote work folder in which VASP was run


.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/latest/
.. _Workchain: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/workflows.html#work-chains
