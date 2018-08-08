# pylint: disable=attribute-defined-outside-init
"""Structure relaxation workchain."""
import enum
from aiida.common.extendeddicts import AttributeDict
from aiida.work import WorkChain
from aiida.work.workchain import append_, while_
from aiida.orm import WorkflowFactory, Code

from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node, init_input
from aiida_vasp.calcs.workchains.restart import UnexpectedCalculationFailure
from aiida_vasp.calcs.workchains.auxiliary.utils import compare_structures, prepare_process_inputs


class RelaxWorkChain(WorkChain):
    """Structure relaxation workchain."""

    _next_workchain = WorkflowFactory('vasp.verify')

    class AlgoEnum(enum.IntEnum):
        """Encode values for algorithm descriptively in enum."""
        NO_UPDATE = -1
        IONIC_RELAXATION_RMM_DIIS = 1
        IONIC_RELAXATION_CG = 2

    class ModeEnum(enum.IntEnum):
        """Encode values for mode of relaxation descriptively in enum."""
        POS_ONLY = 1
        POS_SHAPE_VOL = 3
        POS_SHAPE = 4
        SHAPE_ONLY = 5

        @classmethod
        def get_from_dof(cls, **kwargs):
            """Get the correct mode of relaxation for the given degrees of freedom."""
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
        super(RelaxWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=(get_data_class('structure'), get_data_class('cif')))
        spec.input('potcar_family', valid_type=get_data_class('str'))
        spec.input('potcar_mapping', valid_type=get_data_class('parameter'))
        spec.input('incar', valid_type=get_data_class('parameter'))
        spec.input('options', valid_type=get_data_class('parameter'))
        spec.input('kpoints', valid_type=get_data_class('array.kpoints'), required=False)
        spec.input('settings', valid_type=get_data_class('parameter'), required=False)
        spec.input('restart.max_iterations', valid_type=get_data_class('int'), required=False)
        spec.input('restart.clean_workdir', valid_type=get_data_class('bool'), required=False)
        spec.input('verify.max_iterations', valid_type=get_data_class('int'), required=False)
        spec.input('verify.clean_workdir', valid_type=get_data_class('bool'), required=False)
        spec.input('relax.incar_add', valid_type=get_data_class('parameter'), required=False)
        spec.input('relax.positions', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', True))
        spec.input('relax.shape', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', False))
        spec.input('relax.volume', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', False))
        spec.input('relax.convergence.on', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', False))
        spec.input('relax.convergence.absolute', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', False))
        spec.input('relax.convergence.max_iterations', valid_type=get_data_class('int'), required=False, default=get_data_node('int', 5))
        spec.input(
            'relax.convergence.shape.lengths', valid_type=get_data_class('float'), required=False,
            default=get_data_node('float', 0.1))  # in cartesian coordinates
        spec.input(
            'relax.convergence.shape.angles', valid_type=get_data_class('float'), required=False,
            default=get_data_node('float', 0.1))  # in degree in the cartesian system
        spec.input(
            'relax.convergence.volume', valid_type=get_data_class('float'), required=False,
            default=get_data_node('float', 0.01))  # in degree in the cartesian system
        spec.input(
            'relax.convergence.positions', valid_type=get_data_class('float'), required=False,
            default=get_data_node('float', 0.01))  # in degree in the cartesian system

        spec.outline(
            cls.initialize,
            while_(cls.run_next_workchains)(
                cls.init_next_workchain,
                cls.run_next_workchain,
                cls.verify_next_workchain,
            ),
            cls.results,
            cls.finalize
        )  # yapf: disable

        spec.output('relaxed_structure', valid_type=get_data_class('structure'))
        spec.output('output_parameters', valid_type=get_data_class('parameter'))

    def _set_ibrion(self, incar):
        if self.inputs.relax.positions.value:
            incar.ibrion = self.AlgoEnum.IONIC_RELAXATION_CG
        else:
            incar.ibrion = self.AlgoEnum.NO_UPDATE

    def _set_isif(self, incar):
        """Set relaxation mode according to the chosen degrees of freedom."""
        incar.isif = self.ModeEnum.get_from_dof(
            positions=self.inputs.relax.positions.value, shape=self.inputs.relax.shape.value, volume=self.inputs.relax.volume.value)

    def _add_overrides(self, incar):
        """Add incar tag overrides, except the ones controlled by other inputs (for provenance)."""
        overrides = AttributeDict({k.lower(): v for k, v in self.inputs.relax.incar_add.get_dict().items()})
        if 'ibrion' in overrides:
            raise ValueError('overriding IBRION not allowed, use relax.xxx inputs to control')
        if 'isif' in overrides:
            raise ValueError('overriding ISIF not allowed, use relax.xxx inputs to control')
        if 'nsw' in overrides:
            if self.inputs.relax.positions.value and overrides.nsw < 1:
                raise ValueError('NSW (num ionic steps) was set to 0 but relaxing positions was requested')
            elif not self.inputs.relax.positions.value and overrides.nsw > 0:
                self.report('NSW (num ionic steps) > 1 but relaxing positions was not requested '
                            '(ionic steps will be performed but ions will not move)')
        incar.update(overrides)

    def _assemble_incar(self):
        """Set incar parameters based on other inputs."""
        incar = AttributeDict()
        self._set_ibrion(incar)
        self._set_isif(incar)
        if 'relax.incar_add' in self.inputs:
            try:
                self._add_overrides(incar)
            except ValueError as err:
                self._fail_compat(exception=err)
        return incar

    def initialize(self):
        """Initialize."""
        self._init_context()
        self._init_inputs()
        self._init_structure()

        return

    def _init_context(self):
        """Store exposed inputs in the context."""
        self.ctx.current_restart_folder = None
        self.ctx.is_converged = False
        self.ctx.iteration = 0
        self.ctx.workchains = []

    def _init_structure(self):
        """Initialize the structure."""
        self.ctx.current_structure = self.inputs.structure

        return

    def _init_inputs(self):
        """Initialize the input."""

        self.ctx.inputs = init_input(self.inputs, exclude='relax')

        return

    def run_next_workchains(self):
        within_max_iterations = bool(self.ctx.iteration < self.inputs.relax.convergence.max_iterations.value)
        return bool(within_max_iterations and not self.ctx.is_converged)

    def init_next_workchain(self):
        """Initialize the next workchain calculation."""

        try:
            self.ctx.inputs
        except AttributeError:
            raise ValueError('No input dictionary was defined in self.ctx.inputs')

        self.ctx.inputs.structure = self.ctx.current_structure

    def run_next_workchain(self):
        """Run the next workchain."""

        inputs = prepare_process_inputs(self.ctx.inputs)

        running = self.submit(self._next_workchain, **inputs)

        self.ctx.iteration += 1

        if hasattr(running, 'pid'):
            self.report('launching {}<{}> '.format(self._next_workchain.__name__, running.pid))
        else:
            # Aiida < 1.0
            self.report('launching {}<{}> '.format(self._next_workchain.__name__, running.pk))

        return self.to_context(workchains=append_(running))

    def verify_next_workchain(self):
        """
        Compare the input and output structures of the most recent relaxation run.

        If volume, shape and ion positions are all within a given threshold, consider the relaxation converged.
        """
        if not self.ctx.workchains:
            self._fail_compat(IndexError('The first iteration finished without ' 'returning a {}'.format(self._next_workchain)))

        workchain = self.ctx.workchains[-1]

        # Adopt exit status from last child workchain (supposed to be successfull)
        next_workchain_exit_status = self.ctx.workchains[-1].exit_status
        if not next_workchain_exit_status:
            self.exit_status = 0
        else:
            self.exit_status = next_workchain_exit_status
            self.report('The child {} returned a non-zero exit status, {} '
                        'inherits exit status {}'.format(self._next_workchain, self.__class__.__name__, next_workchain_exit_status))
            return

        if 'output_structure' not in workchain.out:
            self._fail_compat(
                UnexpectedCalculationFailure(
                    'The VaspBaseWf for the relaxation run did not have an output structure and most likely failed'))

        self.ctx.previous_structure = self.ctx.current_structure
        self.ctx.current_structure = workchain.out.output_structure

        converged = True
        if self.inputs.relax.convergence.on.value:
            comparison = compare_structures(self.ctx.previous_structure, self.ctx.current_structure)
            delta = comparison.absolute if self.inputs.relax.convergence.absolute.value else comparison.relative
            if self.inputs.relax.positions.value:
                converged &= self.check_positions_convergence(delta)
            if self.inputs.relax.volume.value:
                converged &= self.check_volume_convergence(delta)
            if self.inputs.relax.shape.value:
                converged &= self.check_shape_convergence(delta)

        if not converged:
            self.ctx.current_restart_folder = workchain.out.remote_folder
            self.report('VaspBaseWf{} was not converged.'.format(workchain))
        elif self.inputs.relax.convergence.on.value:
            self.report('VaspBaseWf{} was converged, finishing.'.format(workchain))
        else:
            self.report('Convergence checking is off')

        self.ctx.is_converged = converged

    def check_shape_convergence(self, delta):
        """Check the difference in cell shape before / after the last iteratio against a tolerance."""
        lengths_converged = bool(delta.cell_lengths.max() <= self.inputs.relax.convergence.shape.lengths.value)
        if not lengths_converged:
            self.report('cell lengths changed by max {}, tolerance is {}'.format(delta.cell_lengths.max(),
                                                                                 self.inputs.relax.convergence.shape.lengths.value))

        angles_converged = bool(delta.cell_angles.max() <= self.inputs.relax.convergence.shape.angles.value)
        if not angles_converged:
            self.report('cell angles changed by max {}, tolerance is {}'.format(delta.cell_angles.max(),
                                                                                self.inputs.relax.convergence.shape.angles.value))

        return bool(lengths_converged and angles_converged)

    def check_volume_convergence(self, delta):
        volume_converged = bool(delta.volume <= self.inputs.relax.convergence.volume.value)
        if not volume_converged:
            self.report('cell volume changed by {}, tolerance is {}'.format(delta.volume, self.inputs.relax.convergence.volume.value))
        return volume_converged

    def check_positions_convergence(self, delta):
        positions_converged = bool(delta.pos_lengths.max() <= self.inputs.relax.convergence.positions.value)
        if not positions_converged:
            self.report('max site position change is {}, tolerance is {}'.format(delta.pos_lengths.max(),
                                                                                 self.inputs.relax.convergence.positions.value))
        return positions_converged

    def results(self):
        """Gather results from the last relaxation iteration."""
        if not self.exit_status:
            last_workchain = self.ctx.workchains[-1]
            relaxed_structure = last_workchain.out.output_structure
            output_parameters = last_workchain.out.output_parameters
            self.out('relaxed_structure', relaxed_structure)
            self.out('output_parameters', output_parameters)
            self.report('workchain completed successfully after {} convergence iterations'.format(self.ctx.iteration))

    def finalize(self):
        """Finalize the workchain."""
        return self.exit_status

    def _fail_compat(self, *args, **kwargs):
        if hasattr(self, 'fail'):
            self.fail(*args, **kwargs)  # pylint: disable=no-member
        else:
            msg = '{}'.format(kwargs['exception'])
            self.abort_nowait(msg)  # pylint: disable=no-member
