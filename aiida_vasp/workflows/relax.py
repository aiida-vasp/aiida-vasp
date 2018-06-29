"""Structure relaxation workchain for VASP."""
import enum
import numpy
from aiida.common.extendeddicts import AttributeDict
from aiida.work import WorkChain
from aiida.work.workchain import append_, while_, ToContext

from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.workflows.base import VaspBaseWf
from aiida_vasp.workflows.restart import prepare_process_inputs, UnexpectedCalculationFailure


def compare_structures(structure_a, structure_b):
    """Compare two StructreData objects A, B and return a delta (A - B) of the relevant properties."""
    delta = AttributeDict()
    volume_a = structure_a.get_cell_volume()
    volume_b = structure_b.get_cell_volume()
    delta.volume = numpy.absolute(volume_a - volume_b)

    pos_a = numpy.array([site.position for site in structure_a.sites])
    pos_b = numpy.array([site.position for site in structure_b.sites])
    delta.pos = pos_a - pos_b

    site_vectors = [delta.pos[i, :] for i in range(delta.pos.shape[0])]
    delta.pos_lengths = numpy.array([numpy.linalg.norm(vector) for vector in site_vectors])

    delta.cell_lengths = numpy.absolute(numpy.array(structure_a.cell_lengths) - numpy.array(structure_b.cell_lengths))

    delta.cell_angles = numpy.absolute(numpy.array(structure_a.cell_angles) - numpy.array(structure_b.cell_angles))

    return delta


class VaspRelaxWf(WorkChain):
    """Structure relaxation workchain for VASP."""

    class IbrionEnum(enum.IntEnum):
        """Encode IBRION values descriptively in enum."""
        NO_UPDATE = -1
        IONIC_RELAXATION_RMM_DIIS = 1
        IONIC_RELAXATION_CG = 2

    class IsifEnum(enum.IntEnum):
        """Encode ISIF values descriptively in enum."""
        POS_ONLY = 1
        POS_SHAPE_VOL = 3
        POS_SHAPE = 4
        SHAPE_ONLY = 5

        @classmethod
        def get_from_dof(cls, **kwargs):
            """Get the correct ISIF value for relaxation for the given degrees of freedom."""
            dof = tuple(kwargs[i] for i in ['positions', 'shape', 'volume'])
            value_from_dof = {
                (True, False, False): cls.POS_ONLY,
                (True, True, True): cls.POS_SHAPE_VOL,
                (True, True, False): cls.POS_SHAPE,
                (False, True, False): cls.SHAPE_ONLY
            }
            return value_from_dof[dof]

    @classmethod
    def define(cls, spec):
        super(VaspRelaxWf, cls).define(spec)
        spec.expose_inputs(VaspBaseWf, include=['code', 'structure', 'potcar_family', 'potcar_mapping', 'options'])
        spec.input('kpoints.mesh', valid_type=get_data_class('array.kpoints'), required=False)
        spec.input('kpoints.distance', valid_type=get_data_class('float'), required=False)
        spec.input('incar_add', valid_type=get_data_class('parameter'), required=False)
        spec.input('relax.positions', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', True))
        spec.input('relax.shape', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', False))
        spec.input('relax.volume', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', False))
        spec.input('convergence.on', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', False))
        spec.input('convergence.max_iterations', valid_type=get_data_class('int'), required=False, default=get_data_node('int', 5))
        spec.input(
            'convergence.shape.lengths', valid_type=get_data_class('float'), required=False,
            default=get_data_node('float', 0.1))  # in cartesian coordinates
        spec.input(
            'convergence.shape.angles', valid_type=get_data_class('float'), required=False,
            default=get_data_node('float', 0.1))  # in degree in the cartesian system
        spec.input(
            'convergence.volume', valid_type=get_data_class('float'), required=False,
            default=get_data_node('float', 0.01))  # in degree in the cartesian system
        spec.input(
            'convergence.positions', valid_type=get_data_class('float'), required=False,
            default=get_data_node('float', 0.01))  # in degree in the cartesian system
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

    def _set_ibrion(self, incar):
        if self.inputs.relax.positions.value:
            incar.ibrion = self.IbrionEnum.IONIC_RELAXATION_CG
        else:
            incar.ibrion = self.IbrionEnum.NO_UPDATE

    def _set_isif(self, incar):
        """Set ISIF value according to the chosen degrees of freedom."""
        incar.isif = self.IsifEnum.get_from_dof(
            positions=self.inputs.relax.positions.value, shape=self.inputs.relax.shape.value, volume=self.inputs.relax.volume.value)

    def _add_overrides(self, incar):
        """Add incar tag overrides, except the ones controlled by other inputs (for provenance)."""
        overrides = AttributeDict({k.lower(): v for k, v in self.inputs.incar_add.get_dict().items()})
        if 'ibrion' in overrides:
            raise ValueError('overriding IBRION not allowed, use relax.xxx inputs to control')
        if 'isif' in overrides:
            raise ValueError('overriding ISIF not allowed, use relax.xxx inputs to control')
        if 'nsw' in overrides:
            if self.inputs.relax.positions.value and overrides.nsw < 1:
                raise ValueError('NSW (num ionic steps) was set to 0 but relaxing positions was requested')
        incar.update(overrides)

    def _assemble_incar(self):
        """Set incar parameters based on other inputs."""
        incar = AttributeDict()
        self._set_ibrion(incar)
        self._set_isif(incar)
        if 'incar_add' in self.inputs:
            try:
                self._add_overrides(incar)
            except ValueError as err:
                self._fail_compat(exception=err)
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
        self.ctx.inputs.incar = self._assemble_incar()
        self.ctx.inputs.code = self.inputs.code
        self.ctx.inputs.structure = self.inputs.structure
        self.ctx.inputs.potcar_family = self.inputs.potcar_family
        self.ctx.inputs.potcar_mapping = self.inputs.potcar_mapping
        self.ctx.inputs.options = self.inputs.options
        if 'max_iterations' in self.inputs.restart:
            self.ctx.inputs.max_iterations = self.inputs.restart.max_iterations

    def validate_inputs(self):
        self.ctx.inputs.kpoints = self._clean_kpoints()

    def should_relax(self):
        within_max_iterations = bool(self.ctx.iteration < self.inputs.convergence.max_iterations.value)
        return bool(within_max_iterations and not self.ctx.is_converged)

    def run_relax(self):
        """Run the BaseVaspWf for the relaxation."""
        if self.ctx.current_restart_folder:
            self.ctx.inputs.restart_folder = self.ctx.current_restart_folder
        self.ctx.inputs.structure = self.ctx.current_structure
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
            self.report('VaspBaseWf{} was converged, finishing.'.format(workchain))
        else:
            self.report('Convergence checking is off')

        self.ctx.is_converged = converged

    def check_shape_convergence(self, delta):
        """Check the difference in cell shape before / after the last iteratio against a tolerance."""
        lengths_converged = bool(delta.cell_lengths.max() <= self.inputs.convergence.shape.lengths.value)
        if not lengths_converged:
            self.report('cell lengths changed by max {}, tolerance is {}'.format(delta.cell_lengths.max(),
                                                                                 self.inputs.convergence.shape.lengths.value))

        angles_converged = bool(delta.cell_angles.max() <= self.inputs.convergence.shape.angles.value)
        if not angles_converged:
            self.report('cell angles changed by max {}, tolerance is {}'.format(delta.cell_angles.max(),
                                                                                self.inputs.convergence.shape.angles.value))

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
            self.fail(*args, **kwargs)  # pylint: disable=no-member
        else:
            msg = '{}'.format(kwargs['exception'])
            self.abort_nowait(msg)  # pylint: disable=no-member
