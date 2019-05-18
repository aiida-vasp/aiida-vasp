# -*- coding: utf-8 -*-
"""Workchain subclass with utilities for restarting, taken from AiiDA-QuantumEspresso."""
from collections import namedtuple
from functools import wraps

from aiida.engine.processes.process import ProcessState
from aiida.common.exceptions import AiidaException
from aiida.common.links import LinkType
from aiida.engine import CalcJob
from aiida.engine import WorkChain, append_
from aiida_vasp.utils.workchains import prepare_process_inputs


class BaseRestartWorkChain(WorkChain):
    """
    Base restart workchain

    This workchain serves as the starting point for more complex workchains that will be designed to
    run a calculation that might need multiple restarts to come to a successful end. These restarts
    may be necessary because a single calculation run is not sufficient to achieve a fully converged
    result, or certain errors maybe encountered which are recoverable.

    This workchain implements the most basic functionality to achieve this goal. It will launch
    calculations, restarting until it is completed successfully or the maximum number of iterations
    is reached. It can recover from simple submission failures and will handle any other errors through
    any of the error handlers that have been provided.

    The idea is to subclass this workchain and leverage the generic error handling that is implemented
    in the few outline methods. Part of the suggested outline would look something like the following::

        cls.init_context
        while_(cls.run_calculations)(
            cls.run_calculation,
            cls.verify_calculation,
        )

    Each of these methods can of course be overriden but they should be general enough to fit most
    calculation cycles. The run_calculation method will take the inputs for the calculation process
    from the context under the key 'inputs'. The user should therefore make sure that before the
    run_calculation method is called, that the to be used inputs are stored under self.ctx.inputs.
    One can update the inputs based on the results from a prior calculation by calling an outline
    method just before the run_calculation step, for example::

        cls.init_context
        while_(cls.run_calculations)(
            cls.prepare_calculation,
            cls.run_calculation,
            cls.verify_calculation,
        )

    Where in the prepare_calculation method, the inputs dictionary at self.ctx.inputs is updated
    before the next calculation will be run with those inputs.
    """
    _verbose = False
    _calculation = None
    _error_handler_entry_point = None
    _expected_calculation_states = [ProcessState.FINISHED, ProcessState.EXCEPTED, ProcessState.KILLED]

    def __init__(self, *args, **kwargs):
        super(BaseRestartWorkChain, self).__init__(*args, **kwargs)
        # Set exit status to None, as this triggers an error if we do not enter
        # a method that detects and sets a particular status
        self.exit_status = None
        if self._calculation is None or not issubclass(self._calculation, CalcJob):
            raise ValueError('no valid CalcJob class defined for _calculation attribute')

        self._error_handlers = []

        return

    @classmethod
    def define(cls, spec):
        super(BaseRestartWorkChain, cls).define(spec)
        spec.exit_code(0, 'NO_ERROR',
            message='the sun is shining')       
        spec.exit_code(100, 'ERROR_ITERATION_RETURNED_NO_CALCULATION',
            message='the run_calculation step did not successfully add a calculation node to the context')
        spec.exit_code(101, 'ERROR_MAXIMUM_ITERATIONS_EXCEEDED',
            message='the maximum number of iterations was exceeded')
        spec.exit_code(102, 'ERROR_UNEXPECTED_CALCULATION_STATE',
            message='the calculation finished with an unexpected calculation state')
        spec.exit_code(103, 'ERROR_SECOND_CONSECUTIVE_SUBMISSION_FAILURE',
            message='the calculation failed to submit, twice in a row')
        spec.exit_code(104, 'ERROR_SECOND_CONSECUTIVE_UNHANDLED_FAILURE',
            message='the calculation failed for an unknown reason, twice in a row')
        spec.exit_code(197, 'ERROR_NO_SOLUTION_FROM_ERROR_HANDLERS',
            message='unknown error detected in the restart workchain')
        spec.exit_code(198, 'ERROR_NO_ERROR_HANDLERS',
            message='unknown error detected in the restart workchain')
        spec.exit_code(199, 'ERROR_UNKNOWN',
            message='unknown error detected in the restart workchain')
        
    def init_context(self):
        """Initialize context variables that are used during the logical flow of the BaseRestartWorkChain."""
        self.ctx.unexpected_failure = False
        self.ctx.submission_failure = False
        self.ctx.restart_calc = None
        self.ctx.is_finished = False
        self.ctx.iteration = 0

        return

    def run_calculations(self):
        """
        Return True whether a new calculation should be run.

        This is the case as long as the last calculation has not finished successfully and the maximum number of restarts
        has not yet been exceeded.
        """
        return not self.ctx.is_finished and self.ctx.iteration < self.inputs.max_iterations.value

    def run_calculation(self):
        """Run a new calculation, taking the input dictionary from the context at self.ctx.inputs."""

        self.ctx.iteration += 1

        try:
            unwrapped_inputs = self.ctx.inputs
        except AttributeError:
            raise ValueError('no calculation input dictionary was defined in self.ctx.inputs')

        inputs = prepare_process_inputs(unwrapped_inputs)
        running = self.submit(self._calculation, **inputs)

        self.report('launching {}<{}> iteration #{}'.format(self._calculation.__name__, running.pk, self.ctx.iteration))

        return self.to_context(calculations=append_(running))

    def verify_calculation(self):
        """
        Analyse the results of the last calculation.

        In particular, check whether it finished successfully or if not troubleshoot the cause and handle the errors,
        or abort if unrecoverable error was found.
        Note that this workchain only check system level issues and does not handle errors or warnings from the
        calculation itself. That is handled in its parent.
        """

        try:
            calculation = self.ctx.calculations[-1]
        except IndexError:
            self.report = 'The first iteration finished without returning a {}'.format(self._calculation.__name__)
            self.exit_status = self.exit_codes.ERROR_ITERATION_RETURNED_NO_CALCULATION
            return self.exit_status
            
        # Set default exit status to an unknown failure
        self.exit_status = self.exit_codes.ERROR_UNKNOWN

        # Done: successful completion of last calculation
        if calculation.is_finished_ok:
            self._handle_succesfull(calculation)

        # Abort: exceeded maximum number of retries
        elif self.ctx.iteration > self.inputs.max_iterations.value:
            self._handle_max_iterations(calculation)

        # Abort: unexpected state of last calculation
        elif calculation.process_state not in self._expected_calculation_states:
            self._handle_unexpected(calculation)

        # Abort: killed
        elif calculation.process_state == ProcessState.KILLED:
            self._handle_killed(calculation)

        # Abort: excepted
        elif calculation.process_state == ProcessState.EXCEPTED:
            self._handle_excepted(calculation)
        
        # Retry or abort: calculation finished or failed (but is not ok)
        elif calculation.process_state == ProcessState.FINISHED:
            self._handle_other(calculation)

        if self.exit_status:
            return self.exit_status

        return self.exit_codes.NO_ERROR

    def results(self):
        """Attach the outputs specified in the output specification from the last completed calculation."""
        if not self.exit_status:
            self.report('{}<{}> completed after {} iterations'.format(self.__class__.__name__, self.pid, self.ctx.iteration))
            
            for name, port in self.spec().outputs.iteritems():
                if port.required and name not in self.ctx.restart_calc.outputs:
                    self.report('the spec specifies the output {} as required '
                                'but was not an output of {}<{}>'.format(name, self._calculation.__name__, self.ctx.restart_calc.pk))

                if name in self.ctx.restart_calc.outputs:
                    node = self.ctx.restart_calc.outputs[name]
                    self.out(name, self.ctx.restart_calc.outputs[name])
                    if self._verbose:
                        self.report("attaching the node {}<{}> as '{}'".format(node.__class__.__name__, node.pk, name))
        else:
            self.report('The calculation {}<{}> returned a non-zero exit status. '
                        'The exit status of {} is thus set to {}'.format(self._calculation, self.pid, self.__class__.__name__,
                                                                         self.exit_status))
        return

    def on_terminated(self):
        """Clean remote folders of the calculations called in the workchain if the clean_workdir input is True."""

        super(BaseRestartWorkChain, self).on_terminated()  # pylint: disable=no-member
        # Do not clean if we do not want to or the calculation failed
        if self.exit_status or self.inputs.clean_workdir.value is False:
            self.report('not cleaning the remote folders')
            return

        cleaned_calcs = []

        for calculation in self.ctx.calculations.outputs(link_type=LinkType.CALL):
            try:
                calculation.outputs.remote_folder._clean()  # pylint: disable=protected-access
                cleaned_calcs.append(calculation.pk)
            except BaseException:
                pass

        if cleaned_calcs:
            self.report('cleaned remote folders of calculations: {}'.format(' '.join(map(str, cleaned_calcs))))

    def _handle_succesfull(self, calculation):
        """Handle the case when the calculaton was successfull."""
        self.report('{}<{}> completed successfully'.format(self._calculation.__name__, calculation.pk))
        self.ctx.restart_calc = calculation
        self.ctx.is_finished = True
        self.exit_status = self.exit_codes.NO_ERROR

        return

    def _handle_max_iterations(self, calculation):
        """Handle the case when the maximum number of iterations are reached."""
        self.report(
            'reached the maximumm number of iterations {}: last ran calculation was {}<{}>'.format(
                self.inputs.max_iterations.value, self._calculation.__name__, calculation.pk))
        self.exit_status = self.exit_codes.ERROR_MAXIMUM_ITERATIONS_EXCEEDED

        return

    def _handle_unexpected(self, calculation):
        """Handle the case when an unexpected state is detected."""
        self.report('unexpected state ({}) of {}<{}>'.format(
            calculation.process_state, self._calculation.__name__, calculation.pk))
        self.exit_status = self.exit_codes.ERROR_UNEXPECTED_CALCULATION_STATE

        return

    def _handle_killed(self, calculation):
        """Handle the case when a killed state is detected, a silent fail."""
        self.exit_status = self.exit_codes.NO_ERROR

        return

    def _handle_excepted(self, calculation):
        """Handle the case when an excepted state is detected, a silent fail."""

        return
    
    # pylint: disable=no-self-use,unused-argument
    def _handle_calculation_sanity_checks(self, calculation):
        """
        Perform additional sanity checks on successfully completed calculation.

        Calculations that run successfully may still have problems that can only be determined when inspecting
        the output. The same problems may also be the hidden root of a calculation failure. For that reason,
        after verifying that the calculation ran, regardless of its calculation state, we perform some sanity
        checks. However, make sure that this workchain is general. Calculation plugin specific sanity checks
        should go into the childs..
        """
        return

    def _handle_unexpected_failure(self, calculation):
        """
        Handle unexpected failure of a calculation.

        The calculation has failed for an unknown reason and could not be handled. If the unexpected_failure
        flag is true, this is the second consecutive unexpected failure and we abort the workchain.
        Otherwise we restart once more.
        """
        if self.ctx.unexpected_failure:
            self.report('failure of {}<{}> could not be handled for the second consecutive time'.format(self._calculation.__name__, calculation.pk))
            self.exit_status = self.exit_codes.ERROR_SECOND_CONSECUTIVE_UNHANDLED_FAILURE
        else:
            self.report('failure of {}<{}> could not be handled, trying to restart'.format(self._calculation.__name__, calculation.pk))
            self.exit_status = self.exit_codes.NO_ERROR

        return

    def _handle_calculation_failure(self, calculation):
        """
        Analyze calculation failure and prepare to restart with corrected inputs.

        The calculation has failed so we try to analyze the reason and change the inputs accordingly
        for the next calculation. If the calculation failed, but did so cleanly, we set it as the
        restart_calc, in all other cases we do not replace the restart_calc
        """
        try:
            outputs = calculation.outputs.output_parameters.get_dict()
            _ = outputs['warnings']
            _ = outputs['parser_warnings']
        except (AttributeError, KeyError) as exception:
            raise UnexpectedCalculationFailure(exception)

        is_handled = False
        handler_report = None
        
        # Sort the handlers based on their priority in reverse order
        handlers = sorted(self._error_handlers, key=lambda x: x.priority, reverse=True)
        
        if not handlers:
            self.exit_status = self.exit_codes.ERROR_NO_ERROR_HANDLERS

        for handler in handlers:
            handler_report = handler.method(self, calculation)

            # If at least one error is handled, we consider the calculation
            # failure handled
            if handler_report and handler_report.is_handled:
                self.exit_status = NO_ERROR
                self.ctx.restart_calc = calculation
                is_handled = True

            # After certain error handlers, we may want to skip all other error
            # handling
            if handler_report and handler_report.do_break:
                break

        # If none of the executed error handlers reported that they handled an
        # error, the failure reason is unknown
        if not is_handled:
            self.exit_status = self.exit_codes.ERROR_NO_SOLUTION_FROM_ERROR_HANDLERS

        return

    def _handle_other(self, calculation):
        """Handle all other cases, not represented by other _handle_*."""

        # Check output for problems independent on calculation state and that
        # do not trigger parser warnings
        self._handle_calculation_sanity_checks(calculation)

        # Calculation failed, try to inspect error code
        self._handle_calculation_failure(calculation)

        # If there is still errors, we need to just try a new resubmit
        if self.exit_status != self.exit_codes.NO_ERROR:
            self._handle_unexpected_failure(calculation)
            self.ctx.unexpected_failure = True
            self.ctx.restart_calc = calculation

        return

    def finalize(self):
        """Finalize the workchain."""

        return self.exit_status


ErrorHandler = namedtuple('ErrorHandler', 'priority method')
"""
A namedtuple to define an error handler for a :class:`~aiida.work.workchain.WorkChain`.

The priority determines in which order the error handling methods are executed, with
the higher priority being executed first. The method defines an unbound WorkChain method
that takes an instance of a :class:`~aiida.orm.implementation.general.calculation.job.AbstractCalcJob`
as its sole argument. If the condition of the error handler is met, it should return an :class:`.ErrorHandlerReport`.

:param priority: integer denoting the error handlers priority
:param method: the workchain class method
"""

ErrorHandlerReport = namedtuple('ErrorHandlerReport', 'is_handled do_break')
"""
A namedtuple to define an error handler report for a :class:`~aiida.work.workchain.WorkChain`.

This namedtuple should be returned by an error handling method of a workchain instance if
the condition of the error handling was met by the failure mode of the calculation.
If the error was appriopriately handled, the 'is_handled' field should be set to `True`,
and `False` otherwise. If no further error handling should be performed after this method
the 'do_break' field should be set to `True`

:param is_handled: boolean, set to `True` when an error was handled
:param do_break: boolean, set to `True` if no further error handling should be performed
"""


# pylint: disable=protected-access
def register_error_handler(cls, priority):
    """
    Turn any decorated function in an error handler for workchain that inherits from the :class:`.BaseRestartWorkChain`.

    The function expects two arguments, a workchain class and a priortity.
    The decorator will add the function as a class method to the workchain class and add an :class:`.ErrorHandler`
    tuple to the :attr:`.BaseRestartWorkChain._error_handlers` attribute of the workchain. During failed calculation
    handling the :meth:`.inspect_calculation` outline method will call the `_handle_calculation_failure` which will loop
    over all error handler in the :attr:`.BaseRestartWorkChain._error_handlers`, sorted with respect to the priority in
    reverse. If the workchain class defines a :attr:`.BaseRestartWorkChain._verbose` attribute and is set to `True`, a
    report message will be fired when the error handler is executed.

    Requirements on the function signature of error handling functions. The function to which the
    decorator is applied needs to take two arguments:

        * `self`: This is the instance of the workchain itself
        * `calculation`: This is the calculation that failed and needs to be investigated

    The function body should usually consist of a single conditional that checks the calculation if
    the error that it is designed to handle is applicable. Although not required, it is advised that
    the function return an :class:`.ErrorHandlerReport` tuple when its conditional was met. If an error was handled
    it should set `is_handled` to `True`. If no other error handlers should be considered set `do_break` to `True`.

    :param cls: the workchain class to register the error handler with
    :param priority: an integer that defines the order in which registered handlers will be called
        during the handling of a failed calculation. Higher priorities will be handled first
    """

    def error_handler_decorator(handler):
        """Turn a function into an error handler for BaseRestartWorkchain."""

        @wraps(handler)
        def error_handler(self, calculation):
            if hasattr(cls, '_verbose') and cls._verbose:
                self.report('({}){}'.format(priority, handler.__name__))
            return handler(self, calculation)

        setattr(cls, handler.__name__, error_handler)

        if not hasattr(cls, '_error_handlers'):
            cls._error_handlers = []
        cls._error_handlers.append(ErrorHandler(priority, error_handler))

        return error_handler

    return error_handler_decorator
