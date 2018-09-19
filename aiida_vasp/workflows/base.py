"""
Base WorkChain for VASP, Error Handling enriched wrapper around VaspCalculation.

Intended to be reused (launched instead of a VaspCalculation) in all other VASP workflows.
Any validation and / or error handling that applies to *every* VASP run,
should be handled on this level, so that every workflow can profit from it.
Anything related to a subset of use cases must be handled in a subclass.
"""
from aiida.work.workchain import while_
from aiida.work.job_processes import override
from aiida.common.extendeddicts import AttributeDict
from aiida.common.exceptions import NotExistent
from aiida.orm import Code, CalculationFactory

from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node, builder_interface
from aiida_vasp.workflows.restart import BaseRestartWorkChain


def get_vasp_proc():
    vasp_cls = CalculationFactory('vasp.vasp')
    if builder_interface(vasp_cls):
        return vasp_cls
    return vasp_cls.process()


class VaspBaseWf(BaseRestartWorkChain):
    """
    Error handling enriched wrapper around VaspCalculation.

    Deliberately conserves most of the interface (required inputs) of the VaspCalculation class.

    This is intended to be used instead of directly submitting a VaspCalculation, so that future features like
    automatic restarting, error checking etc. can be propagated to higher level workchains automatically by implementing them here.

    Usage::

        from aiida.common.extendeddicts import AttributeDict
        from aiida.work import submit
        basevasp = WorkflowFactory('vasp.base')
        inputs = basevasp.get_inputs_template  ## AiiDA < 1.0.0a1
        inputs = basevasp.get_builder() ## AiiDA >= 1.0.0a1
        inputs = AttributeDict()  ## all versions (no tab-completion)
        ## ... set inputs (documented at aiida.readthedocs.io)
        submit(basevasp, **inputs)

    To see working examples, including generation of input nodes from scratch, please refer to ``examples/run_base_wf.py``
    and ``examples/run_vasp.py``.
    """
    _verbose = True
    _calculation_class = CalculationFactory('vasp.vasp')

    @classmethod
    def define(cls, spec):
        super(VaspBaseWf, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=(get_data_class('structure'), get_data_class('cif')))
        spec.input('kpoints', valid_type=get_data_class('array.kpoints'))
        spec.input('potcar_family', valid_type=get_data_class('str'))
        spec.input('potcar_mapping', valid_type=get_data_class('parameter'))
        spec.input('incar', valid_type=get_data_class('parameter'))
        spec.input('wavecar', valid_type=get_data_class('vasp.wavefun'), required=False)
        spec.input('chgcar', valid_type=get_data_class('vasp.chargedensity'), required=False)
        spec.input('settings', valid_type=get_data_class('parameter'), required=False)
        spec.input('options', valid_type=get_data_class('parameter'))

        spec.outline(
            cls.setup,
            cls.validate_inputs,
            while_(cls.should_run_calculation)(
                cls.prepare_calculation,
                cls.run_calculation,
                cls.inspect_calculation
            ),
            cls.results
        )  ## yapf: disable

        spec.output('output_parameters', valid_type=get_data_class('parameter'))
        spec.output('remote_folder', valid_type=get_data_class('remote'))
        spec.output('retrieved', valid_type=get_data_class('folder'))
        spec.output('output_band', valid_type=get_data_class('array.bands'), required=False)
        spec.output('output_structure', valid_type=get_data_class('structure'), required=False)
        spec.output('output_kpoints', valid_type=get_data_class('array.kpoints'), required=False)

    def prepare_calculation(self):
        """Set the restart folder and set INCAR tags for a restart."""
        if isinstance(self.ctx.restart_calc, self._calculation_class):
            self.ctx.inputs.restart_folder = self.ctx.restart_calc.out.remote_folder
            old_incar = AttributeDict(self.ctx.inputs.incar.get_dict())
            incar = old_incar.copy()
            if 'istart' in incar:
                incar.istart = 1
            if 'icharg' in incar:
                incar.icharg = 1
            if incar != old_incar:
                self.ctx.inputs.incar = get_data_node('parameter', dict=incar)

    def validate_inputs(self):
        """Make sure all the required inputs are there and valid, create input dictionary for calculation."""
        self.ctx.inputs = AttributeDict()
        self.ctx.inputs.code = self.inputs.code
        self.ctx.inputs.structure = self.inputs.structure
        self.ctx.inputs.kpoints = self.inputs.kpoints
        self.ctx.inputs.parameters = self.inputs.incar
        if 'settings' in self.inputs:
            self.ctx.inputs.settings = self.inputs.settings

        ## Verify options
        options = AttributeDict()
        options.computer = self.inputs.code.get_computer()
        options.update(self.inputs.options.get_dict())
        #expected_options = ['computer', 'resources']
        #if options.computer.get_scheduler_type() != 'direct':
        #    expected_options.append('queue_name')
        #for option in expected_options:
        #    if option not in options:
        #        self._fail_compat(exception=ValueError('option {} required but not passed!'.format(option)))
        if builder_interface(CalculationFactory('vasp.vasp')):  ## aiida 1.0.0+ will use this
            self.ctx.inputs.options = options
        else:
            self.ctx.inputs._options = options  ## pylint: disable=protected-access

        ## Verify potcars
        try:
            self.ctx.inputs.potential = get_data_class('vasp.potcar').get_potcars_from_structure(
                structure=self.inputs.structure, family_name=self.inputs.potcar_family.value, mapping=self.inputs.potcar_mapping.get_dict())
        except ValueError as err:
            self._fail_compat(exception=err)
        except NotExistent as err:
            self._fail_compat(exception=err)

    @override
    def on_except(self, exc_info):
        """Handle excepted state."""
        try:
            last_calc = self.ctx.calculations[-1] if self.ctx.calculations else None
            if last_calc is not None:
                self.report('Last calculation: {calc}'.format(calc=repr(last_calc)))
                sched_err = last_calc.out.retrieved.get_file_content('_scheduler-stderr.txt')
                sched_out = last_calc.out.retrieved.get_file_content('_scheduler-stdout.txt')
                self.report('Scheduler output:\n{}'.format(sched_out or ''))
                self.report('Scheduler stderr:\n{}'.format(sched_err or ''))
        except AttributeError:
            self.report('No calculation was found in the context. '
                        'Something really awefull happened. Please inspect messages and act.')

        return super(VaspBaseWf, self).on_except(exc_info)
