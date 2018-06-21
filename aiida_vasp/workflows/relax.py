"""Structure relaxation workchain for VASP."""
from aiida.common.extendeddicts import AttributeDict
from aiida.orm.data.base import Float, Bool
from aiida.work import WorkChain
from aiida.work.workchain import append_, ToContext

from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.utils.vasp.isif import IsifStressFlags, Isif
from aiida_vasp.utils.vasp.ibrion import IbrionFlags, Ibrion
from aiida_vasp.utils.vasp.encut import EncutFlags, Encut
from aiida_vasp.utils.vasp.nsw import Nsw
from aiida_vasp.workflows.base import VaspBaseWf
from aiida_vasp.workflows.restart import prepare_process_inputs, UnexpectedCalculationFailure

RELAXATION_INCAR_TEMPLATE = AttributeDict({
    "ismear": -1,
    "sigma": 0.2,
    "algo": 'NORMAL',
    "ediff": 1E-6,
    "prec": 'Normal',
    "ibrion": 2,
    "nsw": 200,
    "nelmin": 4,
    "ediffg": -0.01,
    "isif": 4
})


def relax_parameters(pos, shape, volume):
    """Set IBRION and ISIF according to the choice of what should be relaxed."""
    if not (pos or volume or shape):
        raise ValueError('relaxation workflow requires at least one of "positions", "cell shape", "cell volume" to be alowed to change.')
    if not pos:
        ibrion = Ibrion(ion_updates=IbrionFlags.NO_UPDATE)
    else:
        ibrion = Ibrion(ion_updates=IbrionFlags.IONIC_RELAXATION_CG)
    isif = Isif(calculate_stress=IsifStressFlags.FULL, vary_positions=pos, vary_cell_shape=shape, vary_cell_volume=volume)
    nsw = Nsw(value=200, ibrion=ibrion)
    return ibrion, isif, nsw


def encut_factor(volume):
    if volume:
        return 1.3
    return 1


def clean_incar_overrides(inputs):
    """Make sure the incar overrides do not clash with parameters necessary for relaxation."""
    overrides = AttributeDict(inputs.incar_add.get_dict())
    if 'ibrion' in overrides:
        raise ValueError('overriding IBRION not allowed, use relax.xxx inputs to control')
    if 'isif' in overrides:
        raise ValueError('overriding ISIF not allowed, use relax.xxx inputs to control')
    if 'nsw' in overrides:
        if inputs.relax.positions.value and overrides.nsw < 1:
            raise ValueError('NSW (num ionic steps) was set to 0 but relaxing positions was requested')

    return overrides


class VaspRelaxWf(WorkChain):
    """Structure relaxation workchain for VASP."""

    @classmethod
    def define(cls, spec):
        super(VaspRelaxWf, cls).define(spec)
        spec.expose_inputs(VaspBaseWf, include=['code', 'structure', 'potcar_family', 'potcar_mapping', 'options'])
        spec.input('kpoints.mesh', valid_type=get_data_class('array.kpoints'), required=False)
        spec.input('kpoints.distance', valid_type=Float, required=False)
        spec.input('incar_add', valid_type=get_data_class('parameter'), required=False)
        spec.input('relax.positions', valid_type=Bool, required=False, default=Bool(True))
        spec.input('relax.shape', valid_type=Bool, required=False, default=Bool(False))
        spec.input('relax.volume', valid_type=Bool, required=False, default=Bool(False))
        spec.expose_inputs(VaspBaseWf, namespace='restart', include=['max_iterations'])

        spec.outline(cls.setup, cls.validate_inputs, cls.run_relax, cls.results)

        spec.output('relaxed_structure', valid_type=get_data_class('structure'))
        spec.output('output_parameters', valid_type=get_data_class('parameter'))

    def _clean_incar(self, base_incar):
        """Update incar parameters based on other inputs."""
        incar = base_incar.copy()
        try:
            ibrion, isif, nsw = relax_parameters(
                pos=self.inputs.relax.positions.value, shape=self.inputs.relax.shape.value, volume=self.inputs.relax.shape.value)
            incar.ibrion = ibrion.value
            incar.isif = isif.value
            incar.nsw = nsw.value
        except ValueError as err:
            self._fail_compat(exception=err)

        encut = Encut(
            strategy=EncutFlags.MAX_ENMAX,
            structure=self.inputs.structure,
            potcar_family=self.inputs.potcar_family.value,
            potcar_mapping=self.inputs.potcar_mapping.get_dict(),
            factor=encut_factor(self.inputs.relax.volume))
        if self.inputs.relax.shape.value:
            incar.encut = encut.value

        if 'incar_add' in self.inputs:
            try:
                incar_overrides = clean_incar_overrides(self.inputs)
                incar.update(incar_overrides)
            except ValueError as err:
                self._fail_compat(exception=err)

        for param in [ibrion, isif, encut]:
            errors, warnings = param.clean(incar)
            for error in errors:
                self.logger.error(error)
            for warning in warnings:
                self.report(warning)
            if errors:
                self._fail_compat(exception=ValueError(errors[0]))

        return incar

    def _clean_kpoints(self):
        """Use or create a kpoints mesh depending on user input."""
        kpoints = None
        if 'mesh' in self.inputs.kpoints:
            kpoints = self.inputs.kpoints.mesh
        elif 'distance' in self.inputs.kpoints:
            kpoints = get_data_node(
                'array.kpoints', cell_from_structure=self.inputs.structure, kpoints_mesh_from_density=self.inputs.kpoints.distance.value)
        else:
            self._fail_compat(exception=ValueError('either kpoints.mesh or kpoints.distance is required'))
        return kpoints

    def setup(self):
        """Store exposed inputs in the context."""
        self.ctx.inputs = AttributeDict()
        self.ctx.inputs.code = self.inputs.code
        self.ctx.inputs.structure = self.inputs.structure
        self.ctx.inputs.potcar_family = self.inputs.potcar_family
        self.ctx.inputs.potcar_mapping = self.inputs.potcar_mapping
        self.ctx.inputs.options = self.inputs.options
        if 'max_iterations' in self.inputs.restart:
            self.ctx.inputs.max_iterations = self.inputs.restart.max_iterations

    def validate_inputs(self):
        self.ctx.inputs.kpoints = self._clean_kpoints()
        self.ctx.inputs.incar = self._clean_incar(RELAXATION_INCAR_TEMPLATE)

    def run_relax(self):
        """Run the BaseVaspWf for the relaxation."""
        inputs = prepare_process_inputs(self.ctx.inputs)
        running = self.submit(VaspBaseWf, **inputs)

        self.report('launching VaspBaseWf'.format(running.pid))

        return ToContext(workchains=append_(running))

    def inspect_relax(self):
        """
        Compare the input and output structures of the most recent relaxation run.

        If volume, shape and ion positions are all within a given threshold, consider the relaxation converged.
        """
        if not self.ctx.workchains:
            self._fail_compat(IndexError('The first iteration finished without returning a VaspBaseWf'))
        workchain = self.ctx.workchains[-1]

        if 'output_structure' not in workchain.out:
            self._fail_compat(
                UnexpectedCalculationFailure(
                    'The VaspBaseWf for the relaxation run did not have an output structure and most likely failed'))
        structure = workchain.out.output_structure

        converged = True
        if self.inputs.meta_convergence.ions.value and self.inputs.relax.ions.value:
            converged &= self.check_positions_convergence()
        if self.inputs.meta_convergence.volume.value and self.inputs.relax.volume.value:
            converged &= self.check_volume_convergence()
        if self.inputs.meta_convergence.shape.value and self.inputs.relax.shape.value:
            converged &= self.check_shape_convergence()

    def results(self):
        last_workchain = self.ctx.workchains[-1]
        relaxed_structure = last_workchain.out.output_structure
        output_parameters = last_workchain.out.output_parameters
        self.out('relaxed_structure', relaxed_structure)
        self.out('output_parameters', output_parameters)
