"""
VASP workchain.

---------------
Contains the VaspWorkChain class definition which uses the BaseRestartWorkChain.
"""
from aiida.engine import while_
from aiida.common.lang import override
#from aiida.engine.job_processes import override
from aiida.common.extendeddicts import AttributeDict
from aiida.common.exceptions import NotExistent
from aiida.plugins import CalculationFactory
from aiida.orm import Code
from aiida.engine.processes.workchains.restart import BaseRestartWorkChain, ProcessHandlerReport, process_handler

from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.utils.workchains import compose_exit_code
#from aiida_vasp.workchains.restart import BaseRestartWorkChain
from aiida_vasp.assistant.parameters import ParametersMassage, inherit_and_merge_parameters


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
    _calculation = CalculationFactory('vasp.vasp')

    @classmethod
    def define(cls, spec):  # pylint: disable=too-many-statements
        super(VaspWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=(get_data_class('structure'), get_data_class('cif')), required=True)
        spec.input('kpoints', valid_type=get_data_class('array.kpoints'), required=True)
        spec.input('potential_family', valid_type=get_data_class('str'), required=True)
        spec.input('potential_mapping', valid_type=get_data_class('dict'), required=True)
        spec.input('parameters', valid_type=get_data_class('dict'), required=True)
        spec.input('options', valid_type=get_data_class('dict'), required=True)
        spec.input('settings', valid_type=get_data_class('dict'), required=False)
        spec.input('wavecar', valid_type=get_data_class('vasp.wavefun'), required=False)
        spec.input('chgcar', valid_type=get_data_class('vasp.chargedensity'), required=False)
        spec.input('restart_folder',
                   valid_type=get_data_class('remote'),
                   required=False,
                   help="""
            The restart folder from a previous workchain run that is going to be used.
            """)
        spec.input('max_iterations',
                   valid_type=get_data_class('int'),
                   required=False,
                   default=lambda: get_data_node('int', 5),
                   help="""
            The maximum number of iterations to perform.
            """)
        spec.input('clean_workdir',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', True),
                   help="""
            If True, clean the work dir upon the completion of a successfull calculation.
            """)
        spec.input('verbose',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
            If True, enable more detailed output during workchain execution.
            """)
        spec.input('dynamics.positions_dof',
                   valid_type=get_data_class('list'),
                   required=False,
                   help="""
            Site dependent flag for selective dynamics when performing relaxation
            """)

        spec.outline(
            cls.setup,
            cls.init_inputs,
            while_(cls.should_run_process)(
                cls.prepare_inputs,
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
            cls.finalize
        )  # yapf: disable

        spec.output('misc', valid_type=get_data_class('dict'))
        spec.output('remote_folder', valid_type=get_data_class('remote'))
        spec.output('retrieved', valid_type=get_data_class('folder'))
        spec.output('structure', valid_type=get_data_class('structure'), required=False)
        spec.output('kpoints', valid_type=get_data_class('array.kpoints'), required=False)
        spec.output('trajectory', valid_type=get_data_class('array.trajectory'), required=False)
        spec.output('chgcar', valid_type=get_data_class('vasp.chargedensity'), required=False)
        spec.output('wavecar', valid_type=get_data_class('vasp.wavefun'), required=False)
        spec.output('bands', valid_type=get_data_class('array.bands'), required=False)
        spec.output('forces', valid_type=get_data_class('array'), required=False)
        spec.output('stress', valid_type=get_data_class('array'), required=False)
        spec.output('dos', valid_type=get_data_class('array'), required=False)
        spec.output('occupancies', valid_type=get_data_class('array'), required=False)
        spec.output('energies', valid_type=get_data_class('array'), required=False)
        spec.output('projectors', valid_type=get_data_class('array'), required=False)
        spec.output('dielectrics', valid_type=get_data_class('array'), required=False)
        spec.output('born_charges', valid_type=get_data_class('array'), required=False)
        spec.output('hessian', valid_type=get_data_class('array'), required=False)
        spec.output('dynmat', valid_type=get_data_class('array'), required=False)
        spec.output('site_magnetization', valid_type=get_data_class('dict'), required=False)
        spec.exit_code(0, 'NO_ERROR', message='the sun is shining')
        spec.exit_code(700, 'ERROR_NO_POTENTIAL_FAMILY_NAME', message='the user did not supply a potential family name')
        spec.exit_code(701, 'ERROR_POTENTIAL_VALUE_ERROR', message='ValueError was returned from get_potcars_from_structure')
        spec.exit_code(702, 'ERROR_POTENTIAL_DO_NOT_EXIST', message='the potential does not exist')
        spec.exit_code(703, 'ERROR_IN_PARAMETER_MASSAGER', message='the exception: {exception} was thrown while massaging the parameters')

        # Copied from the old plugin restart workchain
        spec.exit_code(0, 'NO_ERROR', message='the sun is shining')
        spec.exit_code(400,
                       'ERROR_ITERATION_RETURNED_NO_CALCULATION',
                       message='the run_calculation step did not successfully add a calculation node to the context')
        spec.exit_code(401, 'ERROR_MAXIMUM_ITERATIONS_EXCEEDED', message='the maximum number of iterations was exceeded')
        spec.exit_code(402, 'ERROR_UNEXPECTED_CALCULATION_STATE', message='the calculation finished with an unexpected calculation state')
        spec.exit_code(403, 'ERROR_UNEXPECTED_CALCULATION_FAILURE', message='the calculation experienced and unexpected failure')
        spec.exit_code(404, 'ERROR_SECOND_CONSECUTIVE_SUBMISSION_FAILURE', message='the calculation failed to submit, twice in a row')
        spec.exit_code(405,
                       'ERROR_SECOND_CONSECUTIVE_UNHANDLED_FAILURE',
                       message='the calculation failed for an unknown reason, twice in a row')
        spec.exit_code(300,
                       'ERROR_MISSING_REQUIRED_OUTPUT',
                       message='the calculation is missing at least one required output in the restart workchain')
        spec.exit_code(500, 'ERROR_UNKNOWN', message='unknown error detected in the restart workchain')
        spec.exit_code(501,
                       'ERROR_MANUAL_INTERVENTION_NEEDED',
                       message='Cannot handle the error - inputs are likely need to be revised manually.')
        spec.exit_code(502,
                       'ERROR_CALCULATION_NOT_FINISHED',
                       message='Cannot handle the error - the last calculation did not reach the end of execution.')
        spec.exit_code(503,
                       'ERROR_ELECTRONIC_STRUCTURE_NOT_CONVERGED',
                       message='Cannot handle the error - the last calculation did not reach electronic convergence.')

    def _init_parameters(self):
        """Collect input to the workchain in the converge namespace and put that into the parameters."""

        # At some point we will replace this with possibly input checking using the PortNamespace on
        # a dict parameter type. As such we remove the workchain input parameters as node entities. Much of
        # the following is just a workaround until that is in place in AiiDA core.
        parameters = inherit_and_merge_parameters(self.inputs)

        return parameters

    def prepare_inputs(self):
        """
        Set the restart folder and set parameters tags for a restart.
        NOTE: This step may not be needed once we use handlers
        """
        # Check first if the calling workchain wants a restart in the same folder
        if 'restart_folder' in self.inputs:
            self.ctx.inputs.restart_folder = self.inputs.restart_folder
        # Then check if we the restart workchain wants a restart
        if isinstance(self.ctx.restart_calc, self._calculation):
            self.ctx.inputs.restart_folder = self.ctx.restart_calc.outputs.remote_folder
            old_parameters = AttributeDict(self.ctx.inputs.parameters.get_dict())
            parameters = old_parameters.copy()
            if 'istart' in parameters:
                parameters.istart = 1
            if 'icharg' in parameters:
                parameters.icharg = 1
            if parameters != old_parameters:
                self.ctx.inputs.parameters = get_data_node('dict', dict=parameters)

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
            # Also check if the user supplied additional tags that is not in the supported file.
            settings_dict = self.ctx.inputs.settings.get_dict()
            unsupported_parameters = settings_dict.get('unsupported_parameters', unsupported_parameters)
            skip_parameters_validation = settings_dict.get('skip_parameters_validation', skip_parameters_validation)

        # Perform inputs massage to accommodate generalization in higher lying workchains
        # and set parameters.
        try:
            parameters_massager = ParametersMassage(self.ctx.inputs.parameters,
                                                    unsupported_parameters,
                                                    skip_parameters_validation=skip_parameters_validation)
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
                mapping=self.inputs.potential_mapping.get_dict())
        except ValueError as err:
            return compose_exit_code(self.exit_codes.ERROR_POTENTIAL_VALUE_ERROR.status, str(err))  # pylint: disable=no-member
        except NotExistent as err:
            return compose_exit_code(self.exit_codes.ERROR_POTENTIAL_DO_NOT_EXIST.status, str(err))  # pylint: disable=no-member

        try:
            self._verbose = self.inputs.verbose.value
        except AttributeError:
            pass
        # Set the charge density (chgcar)
        if 'chgcar' in self.inputs:
            self.ctx.inputs.charge_density = self.inputs.chgcar

        # Set the wave functions (wavecar)
        if 'wavecar' in self.inputs:
            self.ctx.inputs.wavefunctions = self.inputs.wavecar

        return self.exit_codes.NO_ERROR  # pylint: disable=no-member

    @override
    def on_except(self, exc_info):
        """Handle excepted state."""
        try:
            last_calc = self.ctx.calculations[-1] if self.ctx.calculations else None
            if last_calc is not None:
                self.report('Last calculation: {calc}'.format(calc=repr(last_calc)))  # pylint: disable=not-callable
                sched_err = last_calc.outputs.retrieved.get_file_content('_scheduler-stderr.txt')
                sched_out = last_calc.outputs.retrieved.get_file_content('_scheduler-stdout.txt')
                self.report('Scheduler output:\n{}'.format(sched_out or ''))  # pylint: disable=not-callable
                self.report('Scheduler stderr:\n{}'.format(sched_err or ''))  # pylint: disable=not-callable
        except AttributeError:
            self.report('No calculation was found in the context. '  # pylint: disable=not-callable
                        'Something really awefull happened. '
                        'Please inspect messages and act.')

        return super(VaspWorkChain, self).on_except(exc_info)

    @process_handler(process_handler=1000)
    def _handle_misc_not_exist(self, node):
        """
        Handle the case where misc output is not avaliable, in which case we cannot do anything for it.
        """
        # Check if the run is converged electronically
        if 'misc' not in node.outputs:
            self.report('Cannot found `misc` outputs.')
            return ProcessHandlerReport(exit_code=self.exit_codes.ERROR_UNKNOWN)  # pylint: disable=no-member
        return None

    @process_handler(priority=0)
    def _handle_calculation_sanity_checks(self, node):  # pylint: disable=no-self-use,unused-argument
        """
        Perform additional sanity checks on successfully completed calculation.

        Calculations that run successfully may still have problems that can only be determined when inspecting
        the output. The same problems may also be the hidden root of a calculation failure. For that reason,
        after verifying that the calculation ran, regardless of its calculation state, we perform some sanity
        checks. However, make sure that this workchain is general. Calculation plugin specific sanity checks
        should go into the childs..

        Note that this is the "final" check, most of the potential problems should have been handled and
        fixed in restart by error handlers.
        """

        misc = node.outputs.misc.get_dict()
        if 'run_status' not in 'misc':
            self.report('`run_status` is not found in misc - cannot verify the integrity of the child calcualtio.')
            return ProcessHandlerReport(exit_code=self.exit_codes.ERROR_UNKNOWN)  # pylint: disable=no-member

        run_status = misc['run_status']

        # Check if the calculation is indeed finished
        if not run_status.get('finished'):
            self.report(f'The child calcualtion {node} did not reach the end of execution.')
            return ProcessHandlerReport(exit_code=self.exit_codes.ERROR_CALCUALTION_NOT_FINISHED)  # pylint: disable=no-member

        # Check that the electronic structure is converged
        if not run_status.get('electronic_converged'):
            self.report(f'The child calcualtion {node} did not have converged electronic structure.')
            return ProcessHandlerReport(exit_code=self.exit_codes.ERROR_ELECTRONIC_STRUCTURE_NOT_CONVERGED)  #pylint: disable=no-member

        return None

    @process_handler(priority=100)
    def _handle_unfinished_calculation(self, node):
        """
        Handled the problem such that the calculation is not finished, e.g. did not reach the
        end of execution.

        If WAVECAR exists, just resubmit the calculation with the restart folder.

        If it is a geometry optimisation, attempt to restart with output structure + WAVECAR.
        """

        misc = node.outputs.misc.get_dict()
        run_status = misc['run_status']

        # No need to act if the process is finished
        if run_status.get('finished'):
            return None

        # Check it is a geometry optimisation
        incar = self.ctx.inputs.parameters
        if incar.get('nsw', -1) > 0:
            if 'structure' not in node.outputs.structure:
                self.report('Performing a geometry optimisation but the output structure is not found.')
                return ProcessHandlerReport(do_break=True, exit_code=self.exit_codes.ERROR_MANUAL_INTERVENTION_NEEDED)  #pylint: disable=no-member
            self.report('Continuing geometry optimisation using the last geometry.')
            self.ctx.inputs.structure = node.outputs.structure
            self._attach_wavecar(node)
        else:
            # Single point calculation - attach the wavecar if possible
            has_wavecar = self._attach_wavecar(node)
            if not has_wavecar:
                # If there is no WAVECAR there is not much we can do.....
                self.report('The calcualtions did not finish and no WAVECAR is avalaible for restart.')
                return ProcessHandlerReport(do_break=True, exit_code=self.exit_codes.ERROR_MANUAL_INTERVENTION_NEEDED)  #pylint: disable=no-member
        return ProcessHandlerReport()

    def _attach_wavecar(self, node):
        """
        Attach WAVECAR for the new calculation if there is one.
        """
        # Check the exists of the WAVECAR
        has_wavecar = False
        for file in node.outputs.remote_folder.listdir_withattributes():
            if file['name'] == 'WAVECAR' and file['attributes'].st_size > 1024 * 16:
                has_wavecar = True
        if has_wavecar:
            self.report('Found WAVECAR in the last calcualtion folder, using the outputs for restart.')
            self.ctx.restart_calc = node
        return has_wavecar
