import numpy as np
from aiida.orm import Bool
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.engine import WorkChain

Dict = DataFactory('dict')
KpointsData = DataFactory("array.kpoints")


class BulkModulusWorkChain(WorkChain):
    """WorkChain to compute bulk modulus using VASP."""

    @classmethod
    def define(cls, spec):
        super(BulkModulusWorkChain, cls).define(spec)
        spec.expose_inputs(WorkflowFactory('vasp.relax'))
        spec.outline(
            cls.run_relax,
            cls.run_two_volumes,
        )

    def run_relax(self):
        self.report("run_relax")
        Workflow = WorkflowFactory('vasp.relax')
        builder = Workflow.get_builder()
        builder.code = self.ctx.inputs.code
        builder.parameters = self.ctx.inputs.parameters
        builder.structure = self.ctx.inputs.structure
        builder.settings = self.ctx.inputs.settings
        builder.potential_family = self.ctx.inputs.potential_family
        builder.potential_mapping = self.ctx.inputs.potential_mapping
        builder.kpoints = self.ctx.inputs.kpoints
        builder.options = self.ctx.inputs.options
        builder.metadata.label = self.ctx.inputs.label + " relax"
        builder.metadata.description = self.ctx.inputs.label + " relax"
        builder.clean_workdir = self.ctx.inputs.clean_workdir
        builder.relax = Bool(True)
        builder.force_cutoff = self.ctx.inputs.force_cutoff
        builder.convergence_on = self.ctx.inputs.convergence_on
        builder.convergence_volume = self.ctx.inputs.convergence_volume
        builder.convergence_max_iterations = self.ctx.inputs.convergence_max_iterations
        builder.relax_parameters = self.ctx.inputs.relax_parameters
        builder.verbose = self.ctx.inputs.verbose
        future = self.submit(builder)
        self.to_context(**{'relax': future})

    def run_two_volumes(self):
        self.report("run_two_volumes")
        for strain, future_name in zip((0.99, 1.01), ('minus', 'plus')):
            Workflow = WorkflowFactory('vasp.verify')
            builder = Workflow.get_builder()
            builder.code = self.ctx.inputs.code
            parameters_dict = self.ctx.inputs.parameters.get_dict()
            relax_parameters_dict = self.ctx.inputs.relax_parameters.get_dict()
            parameters_dict.update(relax_parameters_dict)
            parameters_dict.update({'IBRION': 2,
                                    'ISIF': 4})
            structure = self.ctx['relax'].outputs.structure_relaxed.clone()
            structure.set_cell(np.array(structure.cell) * strain ** (1.0 / 3))
            builder.parameters = Dict(dict=parameters_dict)
            builder.structure = self.ctx.inputs.structure
            builder.settings = self.ctx.inputs.settings
            builder.potential_family = self.ctx.inputs.potential_family
            builder.potential_mapping = self.ctx.inputs.potential_mapping
            builder.kpoints = self.ctx.inputs.kpoints
            builder.options = self.ctx.inputs.options
            label = self.ctx.inputs.label + " %s" % future_name
            builder.metadata.label = label
            builder.metadata.description = label
            builder.clean_workdir = self.ctx.inputs.clean_workdir
            builder.verbose = self.ctx.inputs.verbose
            future = self.submit(builder)
            self.to_context(**{future_name: future})
