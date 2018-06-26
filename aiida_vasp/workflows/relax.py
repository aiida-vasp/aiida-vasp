"""Structure relaxation workchain for VASP."""
import numpy
from aiida.common.extendeddicts import AttributeDict
from aiida.orm.data.base import Float, Bool, Int
from aiida.work import WorkChain
from aiida.work.workchain import append_, while_, ToContext

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


def l2_norm(vector_a):
    """L^2 norm for a line vector"""
    return numpy.sqrt(numpy.dot(vector_a, vector_a.T))


def compare_structures(structure_a, structure_b):
    """Compare two StructreData objects A, B and return a delta (A - B) of the relevant properties."""
    delta = AttributeDict()
    volume_a = structure_a.get_cell_volume()
    volume_b = structure_b.get_cell_volume()
    delta.volume = volume_a - volume_b

    pos_a = numpy.array([site.position for site in structure_a.sites])
    pos_b = numpy.array([site.position for site in structure_b.sites])
    delta.pos = pos_a - pos_b

    site_vectors = [delta.pos[i, :] for i in range(delta.pos.shape[0])]
    delta.pos_lengths = numpy.array([l2_norm(vector) for vector in site_vectors])

    delta.cell = numpy.array(structure_a.cell) - numpy.array(structure_b.cell)
    delta.cell_lengths = numpy.array(structure_a.cell_lengths) - numpy.array(structure_b.cell_lengths)
    delta.cell_angles = numpy.array(structure_a.cell_angles) - numpy.array(structure_b.cell_angles)

    return delta


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
        spec.input('convergence.on', valid_type=Bool, required=False, default=Bool(False))
        spec.input('convergence.max_iterations', valid_type=Int, required=False, default=Int(5))
        spec.input('convergence.shape.lengths', valid_type=Float, required=False, default=Float(0.1))  # in cartesian coordinates
        spec.input('convergence.shape.angles', valid_type=Float, required=False, default=Float(0.1))  # in degree in the cartesian system
        spec.input('convergence.volume', valid_type=Float, required=False, default=Float(0.01))  # in degree in the cartesian system
        spec.input('convergence.positions', valid_type=Float, required=False, default=Float(0.01))  # in degree in the cartesian system
        spec.expose_inputs(VaspBaseWf, namespace='restart', include=['max_iterations'])

        spec.outline(
            cls.setup,
            cls.validate_inputs,
            while_(cls.should_relax)(
                cls.run_relax,
                cls.inspect_relax,
            ),
            cls.results
        )  # yapf: disable

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
        self.ctx.current_structure = self.inputs.structure
        self.ctx.current_restart_folder = None
        self.ctx.is_converged = False
        self.ctx.iteration = 0
        self.ctx.workchains = []

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

    def should_relax(self):
        within_max_iterations = bool(self.ctx.iteration < self.inputs.convergence.max_iterations.value)
        return bool(within_max_iterations and not self.ctx.is_converged)

    def run_relax(self):
        """Run the BaseVaspWf for the relaxation."""
        if self.ctx.current_restart_folder:
            self.ctx.inputs.restart_folder = self.ctx.current_restart_folder
        inputs = prepare_process_inputs(self.ctx.inputs)
        running = self.submit(VaspBaseWf, **inputs)
        self.ctx.iteration += 1

        self.report('launching VaspBaseWf{}'.format(running))

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

        self.ctx.previous_structure = self.ctx.current_structure
        self.ctx.current_structure = workchain.out.output_structure

        converged = True
        if self.inputs.convergence.on.value:
            delta = compare_structures(self.ctx.previous_structure, self.ctx.current_structure)
            if self.inputs.relax.positions.value:
                converged &= self.check_positions_convergence(delta)
            if self.inputs.relax.volume.value:
                converged &= self.check_volume_convergence(delta)
            if self.inputs.relax.shape.value:
                converged &= self.check_shape_convergence(delta)

        if not converged:
            self.ctx.current_restart_folder = workchain.out.remote_folder
            self.report('VaspBaseWf{} was not converged.'.format(workchain))
        elif self.inputs.convergence.on.value:
            self.report('Convergence checking is off')
        else:
            self.report('VaspBaseWf{} was converged, finishing.'.format(workchain))

        self.ctx.is_converged = converged

    def check_shape_convergence(self, delta):
        """Check the difference in cell shape before / after the last iteratio against a tolerance."""
        l2_length_changes = l2_norm(delta.cell_lengths)
        lengths_converged = bool(l2_length_changes <= self.inputs.convergence.shape.lengths.value)
        if not lengths_converged:
            self.report('cell lengths changed by {}, tolerance is {}'.format(l2_length_changes,
                                                                             self.inputs.convergence.shape.lengths.value))

        l2_angle_changes = l2_norm(delta.cell_angles)
        angles_converged = bool(l2_angle_changes <= self.inputs.convergence.shape.angles.value)
        if not angles_converged:
            self.report('cell angles changed by {}, tolerance is {}'.format(l2_angle_changes, self.inputs.convergence.shape.angles.value))

        return bool(lengths_converged and angles_converged)

    def check_volume_convergence(self, delta):
        volume_converged = bool(delta.volume <= self.inputs.convergence.volume.value)
        if not volume_converged:
            self.report('cell volume changed by {}, tolerance is {}'.format(delta.volume, self.inputs.convergence.volume.value))
        return volume_converged

    def check_positions_convergence(self, delta):
        positions_converged = bool(delta.pos_lengths.max() <= self.inputs.convergence.positions.value)
        if not positions_converged:
            self.report('max site position change is {}, tolerance is {}'.format(delta.pos_lengths.max(),
                                                                                 self.inputs.convergence.positions.value))
        return positions_converged

    def results(self):
        """Gather results from the last relaxation iteration."""
        last_workchain = self.ctx.workchains[-1]
        relaxed_structure = last_workchain.out.output_structure
        output_parameters = last_workchain.out.output_parameters
        self.out('relaxed_structure', relaxed_structure)
        self.out('output_parameters', output_parameters)
        self.report('workchain completed successfully after {} convergence iterations'.format(self.ctx.iteration))

    def _fail_compat(self, *args, **kwargs):
        if hasattr(self, 'fail'):
            self.fail(*args, **kwargs)
        else:
            msg = '{}'.format(kwargs['exception'])
            self.abort_nowait(msg)
