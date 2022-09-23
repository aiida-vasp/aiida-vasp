"""
VASP workchain.

---------------
Contains the VaspWorkChain class definition which uses the BaseRestartWorkChain.


Below is a copy of the error handler logic from aiida-core.

    If the process is excepted or killed, the work chain will abort. Otherwise any attached handlers will be called
    in order of their specified priority. If the process was failed and no handler returns a report indicating that
    the error was handled, it is considered an unhandled process failure and the process is relaunched. If this
    happens twice in a row, the work chain is aborted. In the case that at least one handler returned a report the
    following matrix determines the logic that is followed:

        Process  Handler    Handler     Action
        result   report?    exit code
        -----------------------------------------
        Success      yes        == 0     Restart
        Success      yes        != 0     Abort
        Failed       yes        == 0     Restart
        Failed       yes        != 0     Abort

    If no handler returned a report and the process finished successfully, the work chain's work is considered done
    and it will move on to the next step that directly follows the `while` conditional, if there is one defined in
    the outline.

This means that for a handler:

    - No error found - just return None
    - No action taken
        - the error is not recoverable - return with a non-zero error code with do break
        - the error is not recoverable, but other handler may/maynot save it - return with a non-zero code without do break
        - the error is not recoverable, and the workchain should be aborted immediately - non-zero code + do break

    - Action taken
        - the error is fixed in full - return with a zero error code with `do_break=True`
        - the error is not fixed in full - return with a report with `do_break=False` but has a `exit_code`.
          this mean other handlers (with lower priority) must handle it and return a zero error_code.

"""
import math

import numpy as np

from aiida.common.exceptions import InputValidationError, NotExistent
from aiida.common.extendeddicts import AttributeDict
from aiida.common.lang import override
from aiida.engine import while_
from aiida.engine.processes.workchains.restart import (
    BaseRestartWorkChain,
    ProcessHandlerReport,
    WorkChain,
    process_handler,
)
from aiida.orm import CalcJobNode, Code
from aiida.plugins import CalculationFactory

from aiida_vasp.assistant.parameters import ParametersMassage, inherit_and_merge_parameters
from aiida_vasp.calcs.vasp import VaspCalculation
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.utils.workchains import compose_exit_code, site_magnetization_to_magmom

# pylint: disable=no-member


class VaspWorkChain(BaseRestartWorkChain):
    """
    The VASP workchain.

    -------------------
    Error handling enriched wrapper around VaspCalculation.

    Deliberately conserves most of the interface (required inputs) of the VaspCalculation class, but
    makes it possible for a user to interact with a workchain and not a calculation.

    This is intended to be used instead of directly submitting a VaspCalculation,
    so that future features like
    automatic restarting, error checking etc. can be propagated to higher level workchains
    automatically by implementing them here.

    Handlers are implemented to try fix common problems and improves the robustness.
    Individual handlers can be enabled/disabled by setting the ``handler_overrides`` input port.
    Additional settings may be passed under the "settings" input, which is also forwarded to the
    calculations. The available options are:

    - ``USE_WAVECAR_FOR_RESTART`` wether calculation restarts should use the WAVECAR. The default is ``True``.

    Usage::

        from aiida.common.extendeddicts import AttributeDict
        from aiida.work import submit
        basevasp = WorkflowFactory('vasp.vasp')
        inputs = basevasp.get_builder()
        inputs = AttributeDict()
        ## ... set inputs
        submit(basevasp, **inputs)

    To see a working example, including generation of input nodes from scratch, please
    refer to ``examples/run_vasp_lean.py``.
    """
    _verbose = False
    _process_class = CalculationFactory('vasp.vasp')
    _algo_switching = {
        'normal': ['fast', 'veryfast', 'damped'],
        'fast': ['normal', 'veryfast', 'damped'],
        'veryfast': ['normal', 'fast', 'damped'],
        'damped': ['normal', 'fast', 'veryfast']
    }

    @classmethod
    def define(cls, spec):  # pylint: disable=too-many-statements
        super(VaspWorkChain, cls).define(spec)
        spec.input(
            'code',
            valid_type=Code,
        )
        spec.input(
            'structure',
            valid_type=(get_data_class('core.structure'), get_data_class('core.cif')),
            required=True,
        )
        spec.input(
            'kpoints',
            valid_type=get_data_class('core.array.kpoints'),
            required=True,
        )
        spec.input(
            'potential_family',
            valid_type=get_data_class('core.str'),
            required=True,
        )
        spec.input(
            'potential_mapping',
            valid_type=get_data_class('core.dict'),
            required=True,
        )
        spec.input(
            'parameters',
            valid_type=get_data_class('core.dict'),
            required=True,
        )
        spec.input(
            'options',
            valid_type=get_data_class('core.dict'),
            required=True,
        )
        spec.input(
            'settings',
            valid_type=get_data_class('core.dict'),
            required=False,
        )
        spec.input(
            'wavecar',
            valid_type=get_data_class('vasp.wavefun'),
            required=False,
        )
        spec.input(
            'chgcar',
            valid_type=get_data_class('vasp.chargedensity'),
            required=False,
        )
        spec.input(
            'site_magnetization',
            valid_type=get_data_class('core.dict'),
            required=False,
            help='Site magnetization to be used as MAGMOM',
        )
        spec.input(
            'restart_folder',
            valid_type=get_data_class('core.remote'),
            required=False,
            help="""
            The restart folder from a previous workchain run that is going to be used.
            """,
        )
        spec.input(
            'max_iterations',
            valid_type=get_data_class('core.int'),
            required=False,
            default=lambda: get_data_node('core.int', 5),
            help="""
            The maximum number of iterations to perform.
            """,
        )
        spec.input(
            'clean_workdir',
            valid_type=get_data_class('core.bool'),
            required=False,
            default=lambda: get_data_node('core.bool', True),
            help="""
            If True, clean the work dir upon the completion of a successful calculation.
            """,
        )
        spec.input(
            'verbose',
            valid_type=get_data_class('core.bool'),
            required=False,
            default=lambda: get_data_node('core.bool', False),
            help="""
            If True, enable more detailed output during workchain execution.
            """,
        )
        spec.input(
            'dynamics.positions_dof',
            valid_type=get_data_class('core.list'),
            required=False,
            help="""
            Site dependent flag for selective dynamics when performing relaxation
            """,
        )
        spec.outline(
            cls.setup,
            cls.init_inputs,
            while_(cls.should_run_process)(
                cls.prepare_inputs,
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )  # yapf: disable

        spec.expose_outputs(cls._process_class)

        # Copied from the old plugin restart workchain
        spec.exit_code(
            0,
            'NO_ERROR',
            message='the sun is shining',
        )
        spec.exit_code(
            300,
            'ERROR_MISSING_REQUIRED_OUTPUT',
            message='the calculation is missing at least one required output in the restart workchain',
        )
        spec.exit_code(
            400,
            'ERROR_ITERATION_RETURNED_NO_CALCULATION',
            message='the run_calculation step did not successfully add a calculation node to the context',
        )
        spec.exit_code(
            401,
            'ERROR_MAXIMUM_ITERATIONS_EXCEEDED',
            message='the maximum number of iterations was exceeded',
        )
        spec.exit_code(
            402,
            'ERROR_UNEXPECTED_CALCULATION_STATE',
            message='the calculation finished with an unexpected calculation state',
        )
        spec.exit_code(
            403,
            'ERROR_UNEXPECTED_CALCULATION_FAILURE',
            message='the calculation experienced and unexpected failure',
        )
        spec.exit_code(
            404,
            'ERROR_SECOND_CONSECUTIVE_SUBMISSION_FAILURE',
            message='the calculation failed to submit, twice in a row',
        )
        spec.exit_code(
            405,
            'ERROR_SECOND_CONSECUTIVE_UNHANDLED_FAILURE',
            message='the calculation failed for an unknown reason, twice in a row',
        )
        spec.exit_code(
            500,
            'ERROR_MISSING_CRITICAL_OUTPUT',
            message='Missing critical output for inspecting the status of the calculation.',
        )
        spec.exit_code(
            501,
            'ERROR_OTHER_INTERVENTION_NEEDED',
            message='Cannot handle the error - inputs are likely need to be revised manually. Message: {message}',
        )
        spec.exit_code(
            502,
            'ERROR_CALCULATION_NOT_FINISHED',
            message='Cannot handle the error - the last calculation did not reach the end of execution.',
        )
        spec.exit_code(
            503,
            'ERROR_ELECTRONIC_STRUCTURE_NOT_CONVERGED',
            message='Cannot handle the error - the last calculation did not reach electronic convergence.',
        )
        spec.exit_code(
            504,
            'ERROR_IONIC_RELAXATION_NOT_CONVERGED',
            message='The ionic relaxation is not converged.',
        )
        spec.exit_code(
            505,
            'ERROR_UNCONVERGED_ELECTRONIC_STRUCTURE_IN_RELAX',
            message=
            'At least one of the ionic steps during the relaxation has did not have converged electronic structure.',
        )
        spec.exit_code(
            700,
            'ERROR_NO_POTENTIAL_FAMILY_NAME',
            message='the user did not supply a potential family name',
        )
        spec.exit_code(
            701,
            'ERROR_POTENTIAL_VALUE_ERROR',
            message='ValueError was returned from get_potcars_from_structure',
        )
        spec.exit_code(
            702,
            'ERROR_POTENTIAL_DO_NOT_EXIST',
            message='the potential does not exist',
        )
        spec.exit_code(
            703,
            'ERROR_IN_PARAMETER_MASSAGER',
            message='the exception: {exception} was thrown while massaging the parameters',
        )

    def setup(self):
        super().setup()
        self.ctx.restart_calc = None
        self.ctx.vasp_did_not_execute = False
        self.ctx.last_calc_was_unfinished = False
        self.ctx.use_wavecar = True
        self.ctx.ignore_transient_nelm_breach = False  # Flag for ignoring the NELM breach during the relaxation
        self.ctx.verbose = None
        self.ctx.last_calc_remote_objects = []
        self.ctx.handler = AttributeDict()
        self.ctx.handler.nbands_increase_tries = 0

    def _init_parameters(self):
        """Collect input to the workchain in the converge namespace and put that into the parameters."""

        # At some point we will replace this with possibly input checking using the PortNamespace on
        # a dict parameter type. As such we remove the workchain input parameters as node entities. Much of
        # the following is just a workaround until that is in place in AiiDA core.
        parameters = inherit_and_merge_parameters(self.inputs)

        return parameters

    def verbose_report(self, *args, **kwargs):
        """Send report if self.ctx.verbose is True"""
        if self.ctx.verbose is True:
            self.report(*args, **kwargs)

    def prepare_inputs(self):
        """
        Enforce some settings for the restart folder and set parameters tags for a restart.
        This is called because launching the sub process.

        NOTE: This method should probably be refactored to give more control on what kind
        of restart is needed
        """
        # Check first if the calling workchain wants a restart in the same folder
        if 'restart_folder' in self.inputs:
            self.ctx.inputs.restart_folder = self.inputs.restart_folder

        # Then check if the workchain wants a restart
        if self.ctx.restart_calc and isinstance(self.ctx.restart_calc.process_class, self._process_class):
            self.ctx.inputs.restart_folder = self.ctx.restart_calc.outputs.remote_folder
            old_parameters = AttributeDict(self.ctx.inputs.parameters).copy()
            parameters = old_parameters.copy()
            # Make sure ISTART and ICHARG is set to read the relevant objects - if they exists
            if 'istart' in parameters and 'WAVECAR' in self.ctx.last_calc_remote_objects:
                # Correct in case of istart = 0
                if parameters.istart == 0 and self.ctx.use_wavecar:
                    parameters.istart = 1
            # Not using the WAVECAR - we make sure ISTART is 0
            if not self.ctx.use_wavecar:
                parameters.istart = 0
            if 'icharg' in parameters and 'CHGCAR' in self.ctx.last_calc_remote_objects:
                parameters.icharg = 1
            if parameters != old_parameters:
                self.ctx.inputs.parameters = parameters
                self.report('Enforced ISTART=1 and ICHARG=1 for restarting the calculation.')

        # Reset the list of valid remote objects and the restart calculation
        self.ctx.last_calc_remote_objects = []
        self.ctx.restart_calc = None

    def update_magmom(self, node=None):
        """
        Update magmom from site magnetization information if available

        :param node: Calculation node to be used, defaults to the last launched calculation.
        """
        if self.is_noncollinear:
            self.report('Automatic carrying on magmom for non-collinear magnetism calculation is not implemented.')
            return

        if node is None:
            node = self.ctx.children[-1]

        if 'site_magnetization' in node.outputs:
            try:
                self.ctx.inputs.parameters['magmom'] = site_magnetization_to_magmom(
                    node.outputs.site_magnetization.get_dict()
                )
            except ValueError:
                pass

    def init_inputs(self):  # pylint: disable=too-many-branches, too-many-statements
        """Make sure all the required inputs are there and valid, create input dictionary for calculation."""
        self.ctx.inputs = AttributeDict()
        self.ctx.inputs.parameters = self._init_parameters()
        # Set the code
        self.ctx.inputs.code = self.inputs.code

        # Set the structure (poscar)
        self.ctx.inputs.structure = self.inputs.structure

        # Set the kpoints (kpoints)
        self.ctx.inputs.kpoints = self.inputs.kpoints

        # Set settings
        unsupported_parameters = None
        skip_parameters_validation = False
        if self.inputs.get('settings'):
            self.ctx.inputs.settings = self.inputs.settings
            # Also check if the user supplied additional tags that is not in the supported object.
            settings_dict = self.ctx.inputs.settings.get_dict()
            unsupported_parameters = settings_dict.get('unsupported_parameters', unsupported_parameters)
            skip_parameters_validation = settings_dict.get('skip_parameters_validation', skip_parameters_validation)
            self.ctx.use_wavecar = settings_dict.get('USE_WAVECAR_FOR_RESTART', True)
            # Ensure that the misc - run_status will be available for parsing - otherwise we abort immediately
            # It should be enabled by default.
            parser_settings = settings_dict.get('parser_settings', {})
            if 'add_misc' in parser_settings:
                misc_spec = parser_settings['add_misc']
                if misc_spec is False:
                    raise InputValidationError('The `misc` output must be requested for parsing!')
                if isinstance(misc_spec, list):
                    if 'run_status' not in misc_spec:
                        raise InputValidationError('The quantity `run_status` must be requested for parsing!')

        # Perform inputs massage to accommodate generalization in higher lying workchains
        # and set parameters.
        try:
            parameters_massager = ParametersMassage(
                self.ctx.inputs.parameters,
                unsupported_parameters,
                skip_parameters_validation=skip_parameters_validation
            )
        except Exception as exception:  # pylint: disable=broad-except
            return self.exit_codes.ERROR_IN_PARAMETER_MASSAGER.format(exception=exception)  # pylint: disable=no-member
        try:
            # Only set if they exists
            # Set any INCAR tags
            self.ctx.inputs.parameters = parameters_massager.parameters.incar
            # Set any dynamics input (currently only for selective dynamics, e.g. custom write to POSCAR)
            self.ctx.inputs.dynamics = parameters_massager.parameters.dynamics
            # Here we could set additional override flags, but those are not relevant for this VASP plugin
        except AttributeError:
            pass

        # Set options
        # Options is very special, not storable and should be
        # wrapped in the metadata dictionary, which is also not storable
        # and should contain an entry for options
        if 'options' in self.inputs:
            options = {}
            options.update(self.inputs.options)
            self.ctx.inputs.metadata = {'options': options}
            # Override the parser name if it is supplied by the user.
            parser_name = self.ctx.inputs.metadata['options'].get('parser_name')
            if parser_name:
                self.ctx.inputs.metadata['options']['parser_name'] = parser_name
            # Set MPI to True, unless the user specifies otherwise
            withmpi = self.ctx.inputs.metadata['options'].get('withmpi', True)
            self.ctx.inputs.metadata['options']['withmpi'] = withmpi

        # Carry on site magnetization for initialization
        if 'site_magnetization' in self.inputs and not self.is_noncollinear:
            magmom = site_magnetization_to_magmom(self.inputs.site_magnetization.get_dict())
            assert len(magmom) == len(self.inputs.structure.sites)
            self.ctx.inputs.parameters['magmom'] = magmom

        # Utilise default input/output selections
        self.ctx.inputs.metadata['options']['input_filename'] = 'INCAR'
        self.ctx.inputs.metadata['options']['output_filename'] = 'OUTCAR'

        # Make sure we also bring along any label and description set on the WorkChain to the CalcJob, it if does
        # not exists, set to empty string.
        if 'metadata' in self.inputs:
            label = self.inputs.metadata.get('label', '')
            description = self.inputs.metadata.get('description', '')
            if 'metadata' not in self.ctx.inputs:
                self.ctx.inputs.metadata = {}
            self.ctx.inputs.metadata['label'] = label
            self.ctx.inputs.metadata['description'] = description

        # Verify and set potentials (potcar)
        if not self.inputs.potential_family.value:
            self.report('An empty string for the potential family name was detected.')  # pylint: disable=not-callable
            return self.exit_codes.ERROR_NO_POTENTIAL_FAMILY_NAME  # pylint: disable=no-member
        try:
            self.ctx.inputs.potential = get_data_class('vasp.potcar').get_potcars_from_structure(
                structure=self.inputs.structure,
                family_name=self.inputs.potential_family.value,
                mapping=self.inputs.potential_mapping.get_dict()
            )
        except ValueError as err:
            return compose_exit_code(self.exit_codes.ERROR_POTENTIAL_VALUE_ERROR.status, str(err))  # pylint: disable=no-member
        except NotExistent as err:
            return compose_exit_code(self.exit_codes.ERROR_POTENTIAL_DO_NOT_EXIST.status, str(err))  # pylint: disable=no-member

        # Store verbose parameter in ctx - otherwise it will not work after deserialization
        try:
            self.ctx.verbose = self.inputs.verbose.value
        except AttributeError:
            self.ctx.verbose = self._verbose

        # Set the charge density (chgcar)
        if 'chgcar' in self.inputs:
            self.ctx.inputs.charge_density = self.inputs.chgcar

        # Set the wave functions (wavecar)
        if 'wavecar' in self.inputs:
            self.ctx.inputs.wavefunctions = self.inputs.wavecar

        return self.exit_codes.NO_ERROR  # pylint: disable=no-member

    @property
    def is_noncollinear(self):
        """Check if the calculation is a noncollinear one"""
        return self.ctx.inputs.parameters.get('lnoncollinear') or self.ctx.inputs.parameters.get('lsorbit')

    @override
    def on_except(self, exc_info):
        """Handle excepted state."""
        try:
            last_calc = self.ctx.calculations[-1] if self.ctx.calculations else None
            if last_calc is not None:
                self.report(f'Last calculation: {repr(last_calc)}')  # pylint: disable=not-callable
                sched_err = last_calc.outputs.retrieved.get_file_content('_scheduler-stderr.txt')
                sched_out = last_calc.outputs.retrieved.get_file_content('_scheduler-stdout.txt')
                self.report(f"Scheduler output:\n{sched_out or ''}")  # pylint: disable=not-callable
                self.report(f"Scheduler stderr:\n{sched_err or ''}")  # pylint: disable=not-callable
        except AttributeError:
            self.report(
                'No calculation was found in the context. '  # pylint: disable=not-callable
                'Something really awful happened. '
                'Please inspect messages and act.'
            )

        return super().on_except(exc_info)

    @override
    def on_terminated(self):
        """
        Clean the working directories of all child calculation jobs if `clean_workdir=True` in the inputs and
        the calculation is finished without problem.
        """
        # Directly called the WorkChain method as this method replaces that of the BaseRestartWorkChain
        WorkChain.on_terminated(self)

        if self.inputs.clean_workdir.value is False:  # type: ignore[union-attr]
            self.report('remote folders will not be cleaned')
            return

        if not self.ctx.is_finished:  # type: ignore[union-attr]
            self.report('remote folders will not be cleaned because the workchain finished with error.')
            return

        cleaned_calcs = []

        for called_descendant in self.node.called_descendants:
            if isinstance(called_descendant, CalcJobNode):
                try:
                    called_descendant.outputs.remote_folder._clean()  # pylint: disable=protected-access
                    cleaned_calcs.append(str(called_descendant.pk))
                except (IOError, OSError, KeyError):
                    pass

        if cleaned_calcs:
            self.report(f"cleaned remote folders of calculations: {' '.join(cleaned_calcs)}")

    @process_handler(priority=2000, enabled=False)
    def handler_always_attach_outputs(self, node):
        """
        Handle the case where we attach the outputs even if underlying child calculation ends up
        with some exit status.
        """

        # Only active this error handler at the last iteration
        if node.is_finished_ok or self.ctx.iteration < self.inputs.max_iterations.value:
            return None

        # Attach all outputs from the last workchain
        self.report('At the last iteration - attaching outputs from the last workchain.')
        self.report('WARNING: The attached outputs may contain incorrect results - proceed with caution.')

        # Attach the required outputs defined in the spec
        for name, port in self.spec().outputs.items():

            try:
                output = node.get_outgoing(link_label_filter=name).one().node
            except ValueError:
                if port.required:
                    self.report(f"required output '{name}' was not an output of {self.ctx.process_name}<{node.pk}>")
            else:
                self.out(name, output)

        # Try to get some meaningful exit codes using the sanity check handler
        # Always generate a handler report with do_break so no more handlers will be run and
        # overwrite the error code.
        report = self._calculation_sanity_checks(node)
        if report:
            self.report('Problems during checks of the outputs. The corresponding `exit_code` will be returned.')
            return ProcessHandlerReport(exit_code=report.exit_code, do_break=True)
        return ProcessHandlerReport(exit_code=self.exit_codes.ERROR_MAXIMUM_ITERATION_EXCEEDED, do_break=True)

    @process_handler(priority=1100, exit_codes=VaspCalculation.exit_codes.ERROR_VASP_DID_NOT_EXECUTE)
    def handler_calculation_did_not_run(self, node):
        """Handle the case where the calculation is not performed"""
        if self.ctx.vasp_did_not_execute:
            self.report(f'{node} did not execute, and this is the second time - aborting.')
            return ProcessHandlerReport(
                do_break=True,
                exit_code=self.exit_codes.ERROR_OTHER_INTERVENTION_NEEDED.format(
                    message='VASP executable did not run on the remote computer.'
                )
            )

        self.report(f'{node} did not execute - try again')
        self.ctx.vasp_did_not_execute = True
        return ProcessHandlerReport(do_break=True)

    @process_handler(priority=1000)
    def handler_misc_not_exist(self, node):
        """
        Handle the case where misc output is not available, in which case we cannot do anything for it.
        """
        # Check if the run is converged electronically
        if 'misc' not in node.outputs:
            self.report('Cannot found `misc` outputs - please check the process reports for issues.')
            return ProcessHandlerReport(exit_code=self.exit_codes.ERROR_MISSING_CRITICAL_OUTPUT, do_break=True)  # pylint: disable=no-member
        return None

    @process_handler(priority=910, exit_codes=[VaspCalculation.exit_codes.ERROR_DID_NOT_FINISH])
    def handler_unfinished_calc_ionic(self, node):
        """
        Handled the problem such that the calculation is not finished, e.g. did not reach the
        end of execution.

        If WAVECAR exists, just resubmit the calculation with the restart folder.

        If it is a geometry optimization, attempt to restart with output structure + WAVECAR.
        """

        # Check it is a geometry optimization
        incar = self.ctx.inputs.parameters
        if incar.get('nsw', -1) > 0:
            if 'structure' not in node.outputs:
                self.report('Performing a geometry optimization but the output structure is not found.')
                return ProcessHandlerReport(
                    do_break=True,
                    exit_code=self.exit_codes.ERROR_OTHER_INTERVENTION_NEEDED.format(
                        message='No output structure for restart.'
                    )
                )  #pylint: disable=no-member
            self.report('Continuing geometry optimization using the last geometry.')
            self.ctx.inputs.structure = node.outputs.structure
            self._setup_restart(node)
            self.update_magmom(node)
            return ProcessHandlerReport(do_break=True)
        return None

    @process_handler(priority=799, enabled=False, exit_codes=[VaspCalculation.exit_codes.ERROR_DID_NOT_FINISH])
    def handler_unfinished_calc_ionic_alt(self, node):
        """
        Handled the problem such that the calculation is not finished, e.g. did not reach the
        end of execution.

        If WAVECAR exists, just resubmit the calculation with the restart folder.

        If it is a geometry optimization, attempt to restart with output structure + WAVECAR.
        """

        # Check it is a geometry optimization
        incar = self.ctx.inputs.parameters
        if incar.get('nsw', -1) > 0:
            if 'structure' not in node.outputs:
                self.report('Performing a geometry optimization but the output structure is not found.')
                return ProcessHandlerReport(
                    do_break=True,
                    exit_code=self.exit_codes.ERROR_OTHER_INTERVENTION_NEEDED.format(
                        message='No output structure for restart.'
                    )
                )  #pylint: disable=no-member
            self.report('Continuing geometry optimization using the last geometry.')
            self.ctx.inputs.structure = node.outputs.structure
            self._setup_restart(node)
            self.update_magmom(node)
            return ProcessHandlerReport(do_break=True)
        return None

    @process_handler(priority=798, enabled=False)
    def handler_unfinished_calc_generic_alt(self, node):
        """
        A generic handler for unfinished calculations, we attempt to restart it once.
        """

        # Only act on this specific return code, otherwise we reset the flag
        if node.exit_status != VaspCalculation.exit_codes.ERROR_DID_NOT_FINISH.status:
            self.ctx.last_calc_was_unfinished = False
            return None

        if self.ctx.last_calc_was_unfinished:
            msg = (
                'The last calculation was not completed for the second time, potentially due to insufficient walltime/node failure.'
                'Please revise the resources request and/or input parameters.'
            )
            return ProcessHandlerReport(
                do_break=True, exit_code=self.exit_codes.ERROR_OTHER_INTERVENTION_NEEDED.format(message=msg)
            )  #pylint: disable=no-member
        self.report((
            'The last calculation was not finished - restart using the same set of inputs. '
            'If it was due to transient problem this may fix it, fingers crossed.'
        ))
        self.ctx.last_calc_was_unfinished = True
        return ProcessHandlerReport(do_break=True)

    @process_handler(priority=900)
    def handler_unfinished_calc_generic(self, node):
        """
        A generic handler for unfinished calculations, we attempt to restart it once.
        """

        # Only act on this specific return code, otherwise we reset the flag
        if node.exit_status != VaspCalculation.exit_codes.ERROR_DID_NOT_FINISH.status:
            self.ctx.last_calc_was_unfinished = False
            return None

        if self.ctx.last_calc_was_unfinished:
            msg = (
                'The last calculation was not completed for the second time, potentially due to insufficient walltime/node failure.'
                'Please revise the resources request and/or input parameters.'
            )
            return ProcessHandlerReport(
                do_break=True, exit_code=self.exit_codes.ERROR_OTHER_INTERVENTION_NEEDED.format(message=msg)
            )  #pylint: disable=no-member
        self.report((
            'The last calculation was not finished - restart using the same set of inputs. '
            'If it was due to transient problem this may fix it, fingers crossed.'
        ))
        self.ctx.last_calc_was_unfinished = True
        return ProcessHandlerReport(do_break=True)

    @process_handler(priority=850, enabled=False)
    def ignore_nelm_breach_relax(self, node):
        """
        Not a actual handler but works as a switch to bypass checks for NELM breaches in the middle of an ionic relaxation.
        """
        _ = node
        self.ctx.ignore_transient_nelm_breach = True

    @process_handler(
        priority=800,
        enabled=False,
        exit_codes=[
            VaspCalculation.exit_codes.ERROR_ELECTRONIC_NOT_CONVERGED,
            VaspCalculation.exit_codes.ERROR_IONIC_NOT_CONVERGED, VaspCalculation.exit_codes.ERROR_DID_NOT_FINISH,
            VaspCalculation.exit_codes.ERROR_VASP_CRITICAL_ERROR, VaspCalculation.exit_codes.ERROR_OVERFLOW_IN_XML
        ]
    )
    def handler_electronic_conv_alt(self, node):  # pylint: disable=too-many-return-statements,too-many-branches
        """Handle electronic convergence problem"""
        incar = node.inputs.parameters.get_dict()
        run_status = node.outputs.misc['run_status']
        notifications = node.outputs.misc['notifications']
        nelm = run_status['nelm']
        algo = incar.get('algo', 'normal')

        # In case of ionic convergence problem, we also act if electronic convergence problem has been reported.
        if node.exit_status in [
            VaspCalculation.exit_codes.ERROR_IONIC_NOT_CONVERGED.status, VaspCalculation.exit_codes.ERROR_DID_NOT_FINISH
        ]:
            perform_fix = False
            if run_status['consistent_nelm_breach']:
                self.report(
                    'The NELM limit has been breached in all ionic steps - proceed to take actions for improving convergence.'
                )
                perform_fix = True
            elif run_status['contains_nelm_breach']:
                # Then there are some breaches in the ionic cycles
                if self.ctx.ignore_transient_nelm_breach:
                    self.report(
                        'WARNING: NELM limit breached in some ionic steps but requested to ignore this - no action taken.'
                    )
                    perform_fix = False
                else:
                    self.report(
                        'The NELM limit has been breached in some ionic steps - proceed to take actions for improving convergence.'
                    )
                    perform_fix = True
            if not perform_fix:
                return None

        if node.exit_status == VaspCalculation.exit_codes.ERROR_VASP_CRITICAL_ERROR.status:
            # Make sure we only continue in this handler for a selected set of the critical errors.
            if not any(item in node.exit_message for item in ['EDDRMM', 'EDDDAV', 'The topmost band is occupied']):
                # We have some other critical error not to handle here
                return None

        # Check if we need to add more bands
        for item in notifications:
            if item['name'] == 'bandocc' and self.ctx.handler.nbands_increase_tries < 5:
                try:
                    nbands = run_status['nbands']
                    # Increase nbands with 10%
                    nbands_new = math.ceil(nbands * 1.1)
                    self.report(f'Changing NBANDS from {nbands} to {nbands_new}')
                    self.ctx.handler.nbands_increase_tries += 1
                    incar['nbands'] = nbands_new
                    self.ctx.inputs.parameters.update(incar)
                    return ProcessHandlerReport(do_break=True)
                except KeyError:
                    self.report(
                        'The topmost band is occupied but did not locate the nbands entry in run_status, so no way to do corrections.'
                    )
                    return ProcessHandlerReport(exit_code=self.exit_codes.ERROR_MISSING_CRITICAL_OUTPUT, do_break=True)  # pylint: disable=no-member

        if nelm < 300:
            # Standard NELM might be a bit low, so increase a bit
            incar['nelm'] = 300
            # Here we can just continue from previous run
            self._setup_restart(node)
            self.ctx.inputs.parameters.update(incar)
            self.report(f'Changing NELM from {nelm} to 300.')
            return ProcessHandlerReport(do_break=True)

        # Let us start or continue to switch algorithms and reduce try list
        try:
            if self.ctx.handler.get('remaining_algos', None) is None:
                self.ctx.handler.remaining_algos = self._algo_switching[algo.lower()]
            new_algo = self.ctx.handler.remaining_algos.pop(0)
            incar['algo'] = new_algo
            self.ctx.inputs.parameters.update(incar)
            self.report(f'Changing ALGO from {algo.lower()} to {new_algo}')
            return ProcessHandlerReport(do_break=True)
        except IndexError:
            self.report('No more algorithms to try and we still have not reached electronic convergence.')

        self.report('No additional fixes can be applied to improve the electronic convergence - aborting.')
        return ProcessHandlerReport(
            do_break=True,
            exit_code=self.exit_codes.ERROR_OTHER_INTERVENTION_NEEDED.format(
                message='Cannot apply fix for reaching electronic convergence.'
            )
        )

    @process_handler(
        priority=800,
        exit_codes=[
            VaspCalculation.exit_codes.ERROR_ELECTRONIC_NOT_CONVERGED,
            VaspCalculation.exit_codes.ERROR_IONIC_NOT_CONVERGED,
            VaspCalculation.exit_codes.ERROR_DID_NOT_FINISH,
        ]
    )
    def handler_electronic_conv(self, node):
        """Handle electronic convergence problem"""
        incar = node.inputs.parameters.get_dict()
        run_status = node.outputs.misc['run_status']
        nelm = run_status['nelm']
        algo = incar.get('algo', 'normal')

        # In case of ionic convergence problem, we also act if electronic convergence problem has been reported.
        if node.exit_status in [
            VaspCalculation.exit_codes.ERROR_IONIC_NOT_CONVERGED.status, VaspCalculation.exit_codes.ERROR_DID_NOT_FINISH
        ]:
            perform_fix = False
            if run_status['consistent_nelm_breach']:
                self.report(
                    'The NELM limit has been breached in all ionic steps - proceed to take actions for improving convergence.'
                )
                perform_fix = True
            elif run_status['contains_nelm_breach']:
                # Then there are some breaches in the ionic cycles
                if self.ctx.ignore_transient_nelm_breach:
                    self.report(
                        'WARNING: NELM limit breached in some ionic steps but requested to ignore this - no action taken.'
                    )
                    perform_fix = False
                else:
                    self.report(
                        'The NELM limit has been breached in some ionic steps - proceed to take actions for improving convergence.'
                    )
                    perform_fix = True
            if not perform_fix:
                return None

        if algo.lower() in ('fast', 'veryfast'):
            incar['algo'] = 'normal'
            self._setup_restart(node)
            self.ctx.inputs.parameters.update(incar)
            self.report(f'Setting ALGO=normal from ALGO={algo.lower()}')
            return ProcessHandlerReport(do_break=True)

        # The logic below only works for algo=normal
        if algo.lower() == 'normal':
            # First try - Increase NELM
            if nelm < 150:
                incar['nelm'] = 150
                self._setup_restart(node)
                self.ctx.inputs.parameters.update(incar)
                self.report('Setting NELM to 150')
                return ProcessHandlerReport(do_break=True)
            # Adjust AMIX value if NELM is already high
            amix = incar.get('amix', 0.4)
            amix_steps = [0.2, 0.1, 0.05]
            for amix_target in amix_steps:
                if amix > amix_target:
                    incar['amix'] = amix_target
                    # Increase NELM in the mean time - smaller amplitude requires more cycles but more stable.
                    incar['nelm'] = nelm + 20
                    self._setup_restart(node)
                    self.ctx.inputs.parameters.update(incar)
                    self.report(f"Reducing AMIX to {incar['amix']}")
                    return ProcessHandlerReport(do_break=True)
            # Change to ALGO if options have been exhausted
            incar['algo'] = 'all'
            self.ctx.inputs.parameters.update(incar)
            self._setup_restart(node)
            self.report('Switching to ALGO = ALL')
            return ProcessHandlerReport(do_break=True)
        self.report('No additional fixes can be applied to improve the electronic convergence - aborting.')
        return ProcessHandlerReport(
            do_break=True,
            exit_code=self.exit_codes.ERROR_OTHER_INTERVENTION_NEEDED.format(
                message='Cannot apply fix for reaching electronic convergence.'
            )
        )

    @process_handler(priority=510, exit_codes=[VaspCalculation.exit_codes.ERROR_IONIC_NOT_CONVERGED], enabled=False)
    def handler_ionic_conv_enhanced(self, node):  #pylint: disable=too-many-return-statements, too-many-branches
        """
        Enhanced handling of ionic relaxation problem beyond simple restarts.

        This is only used when the calculation is having difficulties reaching the
        convergence. This handler should be applied before the standard handler which
        breaks the handling cycle.
        """

        if 'structure' not in node.outputs:
            self.report('Performing a geometry optimization but the output structure is not found.')
            return ProcessHandlerReport(
                do_break=True,
                exit_code=self.exit_codes.ERROR_OTHER_INTERVENTION_NEEDED.format(
                    message='No output structure for restarting ionic relaxation.'
                )
            )  #pylint: disable=no-member

        # The simplest solution - resubmit the calculation again
        child_nodes = self.ctx.children
        child_miscs = [node.outputs.misc for node in child_nodes]

        self.update_magmom(node)

        # Enhanced handler only takes place after 3 trials
        if len(child_miscs) < 3:
            return None

        natom = len(self.ctx.inputs.structure.sites)
        # Number of iterations
        ionic_iterations = [misc['run_status']['last_iteration_index'][0] for misc in child_miscs]

        # Output energies
        energies = []
        for misc in child_miscs:
            energies.append(misc.get('total_energies', {}).get('energy_extrapolated'))
        if all([eng is not None for eng in energies[-3:]]):  # pylint: disable=use-a-generator
            de_per_atom = np.diff(energies) / natom
        else:
            return None

        # First check if dE is very small
        if np.all(np.abs(de_per_atom) < 1e-5):
            msg = (
                'The total energy difference between the last two step is smaller than 1e-5 /atom'
                '- please consider to revise the cutoff value of the ionic steps.'
            )
            self.report(msg)
            return ProcessHandlerReport(
                do_break=True, exit_code=self.exit_codes.ERROR_OTHER_INTERVENTION_NEEDED.format(msg)
            )

        # Check if there are very few step performed per launch. Because VASP does not carry over
        # the internal parameters of the optimizer, this can make convergence slower.
        if ionic_iterations[-1] < 5:
            msg = (
                'Less than 5 iterations performed in the last launch - '
                'please consider submitting the jobs with revised resources request.'
            )
            self.report(msg)
            return ProcessHandlerReport(
                do_break=True, exit_code=self.exit_codes.ERROR_OTHER_INTERVENTION_NEEDED.format(msg)
            )

        # Warn about very unusually large number of steps and switch IBRION if needed.
        # Total degrees of freedom
        dof = 3 * natom
        isif = self.ctx.inputs.parameters.get('isif', 2)
        if isif == 3:
            dof += 6

        if sum(ionic_iterations) > dof + 10:
            self.report(f'Unusually large number of iterations performed for the degrees of freedom: {dof}')
            ibrion = self.ctx.inputs.parameters.get('ibrion')
            # In this case alternate between different relaxation algorithms
            if ibrion == 2:
                self.ctx.inputs.parameters['ibrion'] = 1
                self.ctx.inputs.parameters['potim'] = 0.3
                self.ctx.inputs.structure = node.outputs.structure
                self.report('Switching to IBRION=1 from IBRION=2 with POTIM = 0.3')
                return ProcessHandlerReport(do_break=True)
            if ibrion == 1:
                self.ctx.inputs.parameters['ibrion'] = 2
                self.ctx.inputs.parameters['potim'] = 0.1
                self.ctx.inputs.structure = node.outputs.structure
                self.report('Switching to IBRION=2 from IBRION=1 with POTIM = 0.1')
                return ProcessHandlerReport(do_break=True)

        # Check if energies are increasing without significant volume change
        vol_changes = []
        for child in child_nodes:
            inp_vol = child.inputs.structure.get_cell_volume()
            out_vol = child.outputs.structure.get_cell_volume()
            vol_changes.append(out_vol / inp_vol - 1.0)
        vol_changes = np.array(vol_changes)

        vol_tol = 0.03
        if np.all(de_per_atom > 0.0) and np.all(abs(vol_changes[-2:]) < vol_tol):
            msg = 'Energy increasing for the last two iterations - something can be very wrong...'
            self.report(msg)
            return ProcessHandlerReport(
                do_break=True, exit_code=self.exit_codes.ERROR_OTHER_INTERVENTION_NEEDED.format(msg)
            )

        self.report('No fixes can be applied for ionic convergence.')
        return None

    @process_handler(priority=505, exit_codes=[VaspCalculation.exit_codes.ERROR_IONIC_NOT_CONVERGED])
    def handler_ionic_conv(self, node):
        """Handle ionic convergence problem"""
        if 'structure' not in node.outputs:
            self.report('Performing a geometry optimization but the output structure is not found.')
            return ProcessHandlerReport(
                do_break=True,
                exit_code=self.exit_codes.ERROR_OTHER_INTERVENTION_NEEDED.format(
                    message='No output structure for restarting ionic relaxation.'
                )
            )  #pylint: disable=no-member
        # The simplest solution - resubmit the calculation again
        self.report('Continuing geometry optimization using the last geometry.')
        self.ctx.inputs.structure = node.outputs.structure
        self._setup_restart(node)
        self.update_magmom(node)
        return ProcessHandlerReport(do_break=True)

    @process_handler(priority=400, exit_codes=[VaspCalculation.exit_codes.ERROR_VASP_CRITICAL_ERROR])
    def handler_vasp_critical_error(self, node):
        """
        Check if the calculation contain any critical error.
        """
        notification = node.outputs.misc['notifications']
        message = f"Critical error detected in the notifications: {', '.join([item.get('name') for item in notification])}"
        self.report(message + ' - aborting.')
        return ProcessHandlerReport(
            do_break=True, exit_code=self.exit_codes.ERROR_OTHER_INTERVENTION_NEEDED.format(message=message)
        )

    @process_handler(priority=5)
    def check_misc_output(self, node):
        """
        Check if misc output exists.
        """
        misc = node.outputs.misc.get_dict()
        if 'run_status' not in misc:
            self.report('`run_status` is not found in misc - cannot verify the integrity of the child calculation.')
            return ProcessHandlerReport(exit_code=self.exit_codes.ERROR_MISSING_CRITICAL_OUTPUT, do_break=True)  # pylint: disable=no-member
        return None

    @process_handler(priority=4)
    def check_calc_is_finished(self, node):
        """
        Check if the calculation has reached the end of execution.
        """
        misc = node.outputs.misc.get_dict()
        run_status = misc['run_status']
        if not run_status.get('finished'):
            self.report(f'The child calculation {node} did not reach the end of execution.')
            return ProcessHandlerReport(exit_code=self.exit_codes.ERROR_CALCULATION_NOT_FINISHED, do_break=True)  # pylint: disable=no-member
        return None

    @process_handler(priority=3)
    def check_electronic_converged(self, node):
        """
        Check if the calculation has converged electronic structure.
        """
        misc = node.outputs.misc.get_dict()
        run_status = misc['run_status']
        # Check that the electronic structure is converged
        if not run_status.get('electronic_converged'):
            self.report(f'The child calculation {node} does not possess a converged electronic structure.')
            return ProcessHandlerReport(
                exit_code=self.exit_codes.ERROR_ELECTRONIC_STRUCTURE_NOT_CONVERGED, do_break=True
            )  #pylint: disable=no-member
        if run_status.get('contains_nelm_breach'):
            if self.ctx.ignore_transient_nelm_breach:
                self.report(
                    'The calculation contains at least one electronic minimization '
                    'that was truncated. It should thus not be considered converged. '
                    'Upon request from user, this is ignored.'
                )
            else:
                self.report(
                    'The calculation contains at least one electronic minimization '
                    'that is truncated. It should thus not be considered converged. '
                    'Treating the calculation as failed. Please inspect, maybe it is salvageable.'
                )
                return ProcessHandlerReport(
                    exit_code=self.exit_codes.ERROR_UNCONVERGED_ELECTRONIC_STRUCTURE_IN_RELAX, do_break=True
                )  #pylint: disable=no-member

        return None

    @process_handler(priority=2)
    def check_ionic_converged(self, node):
        """
        Check if the calculation has converged ionic structure.
        """

        # Check if we have requested to ignore ionic convergence check at calculation level
        # If so, then this handler should be by-passed
        if 'settings' in node.inputs:
            settings = node.inputs.settings.get_dict()
            if not settings.get('CHECK_IONIC_CONVERGENCE', True):
                return None

        misc = node.outputs.misc.get_dict()
        run_status = misc['run_status']

        # Check that the ionic structure is converged
        if run_status.get('ionic_converged') is False:
            self.report(f'The child calculation {node} did not have converged ionic structure.')
            return ProcessHandlerReport(exit_code=self.exit_codes.ERROR_IONIC_RELAXATION_NOT_CONVERGED, do_break=True)  #pylint: disable=no-member
        return None

    def _calculation_sanity_checks(self, node):  # pylint: disable=unused-argument
        """
        Perform additional sanity checks on successfully completed calculation.
        This method acts invokes the 'check' handlers to check the calculations and abort the workchain if any
        problem is found. This is useful when all of the corresponding error handlers are disabled, and allow
        one to avoid the default behaviour of restarting the calculation one more times regardlessly with unhandled errors.
        """
        checks = [
            self._check_misc_output, self._check_calc_is_finished, self._check_electronic_converged,
            self._check_ionic_converged
        ]

        # Go though the checks one after another, return report if necessary
        last_report = None
        for check in checks:
            report = check(node)
            if report:
                if report.do_break:
                    return report
                last_report = report
        return last_report

    def _update_last_calc_objects(self, node):
        """
        Connect to the remote and find the valid objects in th calculation folder

        Only update if the entry is empty in order to avoid too many connections to the remote.
        """
        if not self.ctx.last_calc_remote_objects:
            self.ctx.last_calc_remote_objects = list_valid_objects_in_remote(node.outputs.remote_folder)
        return self.ctx.last_calc_remote_objects

    def _setup_restart(self, node):
        """
        Check the existence of any restart objects, if any of them eixsts use the last calculation
        for restart.
        """
        self._update_last_calc_objects(node)
        if 'WAVECAR' in self.ctx.last_calc_remote_objects or 'CHGCAR' in self.ctx.last_calc_remote_objects:
            self.ctx.restart_calc = node
            return True
        return False


def list_valid_objects_in_remote(remote, path='.', size_threshold=0) -> list:
    """
    List non-empty objects in the remote folder

    :param remote: The `RemoteFolder` node to be inspected.
    :param path: The relative path.
    :param size_threshold: The size threshold to treat the object as a valide one.

    :returns: A list of valid objects in the directory.
    """
    none_empty = []
    try:
        contents = remote.listdir_withattributes(path)
    except OSError:
        return []

    for obj in contents:
        if obj['attributes'].st_size > size_threshold and not obj['isdir']:
            none_empty.append(obj['name'])
    return none_empty
