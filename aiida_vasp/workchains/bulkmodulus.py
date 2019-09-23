import numpy as np
from aiida.orm import Bool
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.engine import WorkChain, calcfunction
from aiida.common.extendeddicts import AttributeDict

Dict = DataFactory('dict')
Float = DataFactory('float')
KpointsData = DataFactory("array.kpoints")


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
            cls.run_relax,
            cls.run_two_volumes,
            cls.calc_bulk_modulus
        )
        spec.output('bulk_modulus', valid_type=Float)

    def run_relax(self):
        self.ctx.inputs = AttributeDict()
        self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))
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
        builder.metadata.label = self.ctx.inputs.metadata.label + " relax"
        builder.metadata.description = self.ctx.inputs.metadata.label + " relax"
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
            structure = get_strained_structure(
                self.ctx['relax'].outputs.structure_relaxed,
                Float(strain))
            builder.parameters = Dict(dict=parameters_dict)
            builder.structure = structure
            builder.settings = self.ctx.inputs.settings
            builder.potential_family = self.ctx.inputs.potential_family
            builder.potential_mapping = self.ctx.inputs.potential_mapping
            builder.kpoints = self.ctx.inputs.kpoints
            builder.options = self.ctx.inputs.options
            label = self.ctx.inputs.metadata.label + " %s" % future_name
            builder.metadata.label = label
            builder.metadata.description = label
            builder.clean_workdir = self.ctx.inputs.clean_workdir
            builder.verbose = self.ctx.inputs.verbose
            future = self.submit(builder)
            self.to_context(**{future_name: future})

    def calc_bulk_modulus(self):
        self.report("calc_bulk_modulus")
        bulk_modulus = calculate_bulk_modulus(
            self.ctx['minus'].outputs.stress,
            self.ctx['plus'].outputs.stress,
            self.ctx['minus'].inputs.structure,
            self.ctx['plus'].inputs.structure)
        self.out('bulk_modulus', bulk_modulus)
        self.report('finish bulk modulus calculation')
