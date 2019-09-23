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
            cls.initialize,
            cls.run_relax,
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
            structure = get_strained_structure(
                self.ctx['relax'].outputs.structure_relaxed,
                Float(strain))
            builder.structure = structure
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
        self.out('bulk_modulus', bulk_modulus)
        self.report('finish bulk modulus calculation')
