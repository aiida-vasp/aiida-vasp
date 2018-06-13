# -*- coding: utf-8 -*-
"""Workchain subclass with utilities for restarting, taken from AiiDA-QuantumEspresso."""
from collections import namedtuple
from functools import wraps

from aiida.common.datastructures import calc_states
from aiida.common.exceptions import AiidaException
from aiida.common.extendeddicts import AttributeDict
from aiida.common.links import LinkType
from aiida.orm.calculation import JobCalculation
from aiida.orm.data.base import Bool, Int
from aiida.orm.data.parameter import ParameterData
from aiida.work.workchain import WorkChain, ToContext, append_


def prepare_process_inputs(inputs):
    """
    Prepare the inputs dictionary for a calculation process.

    Any remaining bare dictionaries in the inputs dictionary will be wrapped in a ParameterData data node
    except for the '_options' key which should remain a standard dictionary. Another exception are dictionaries
    whose keys are not strings but for example tuples.
    This is the format used by input groups as in for example the explicit pseudo dictionary where the key is
    a tuple of kind to which the UpfData corresponds.
    """
    prepared_inputs = AttributeDict()

    for key, val in inputs.iteritems():
        if key not in ['_options', 'options'] and isinstance(val, dict) and all([isinstance(k, (basestring)) for k in val.keys()]):
            prepared_inputs[key] = ParameterData(dict=val)
        else:
            prepared_inputs[key] = val

    return prepared_inputs


def finished_ok_compat(calc):
    if hasattr(calc, 'has_finished_ok'):
        return calc.has_finished_ok()
    return calc.finished_ok


class UnexpectedCalculationFailure(AiidaException):
    """Raised when a Calculation has failed for an unknown reason."""


class MaxIterationsFailure(AiidaException):
    """Raised when a Calculation has failed after the maximum allowed number of tries."""


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

        cls.setup
        while_(cls.should_run_calculation)(
            cls.run_calculation,
            cls.inspect_calculation,
        )

    Each of these methods can of course be overriden but they should be general enough to fit most
    calculation cycles. The run_calculation method will take the inputs for the calculation process
    from the context under the key 'inputs'. The user should therefore make sure that before the
    run_calculation method is called, that the to be used inputs are stored under self.ctx.inputs.
    One can update the inputs based on the results from a prior calculation by calling an outline
    method just before the run_calculation step, for example::

        cls.setup
        while_(cls.should_run_calculation)(
            cls.prepare_calculation,
            cls.run_calculation,
            cls.inspect_calculation,
        )

    Where in the prepare_calculation method, the inputs dictionary at self.ctx.inputs is updated
    before the next calculation will be run with those inputs.
    """
    _verbose = False
    _calculation_class = None
    _error_handler_entry_point = None
    _expected_calculation_states = [calc_states.FINISHED, calc_states.FAILED, calc_states.SUBMISSIONFAILED]

    def _fail_compat(self, *args, **kwargs):
        if hasattr(self, 'fail'):
            self.fail(*args, **kwargs)
        else:
            msg = '{}'.format(kwargs['exception'])
            self.abort_nowait(msg)  # pylint: disable=no-member

    def __init__(self, *args, **kwargs):
        super(BaseRestartWorkChain, self).__init__(*args, **kwargs)

        if self._calculation_class is None or not issubclass(self._calculation_class, JobCalculation):
            raise ValueError('no valid JobCalculation class defined for _calculation_class attribute')

        ## compatibility v0.12.x <-> v1.0.0-alpha
        self._error_handlers = []
        if not hasattr(super(BaseRestartWorkChain, self), 'on_destroy'):
            self.on_terminated = self.on_destroy
        return

    @classmethod
    def define(cls, spec):
        super(BaseRestartWorkChain, cls).define(spec)
        spec.input(
            'max_iterations',
            valid_type=Int,
            default=Int(5),
            help="""
            the maximum number of iterations the workchain will attempt to get the calculation to finish successfully
            """)
        spec.input(
            'clean_workdir',
            valid_type=Bool,
            default=Bool(False),
            help="""
            when set to True, the work directories of all called calculation will be cleaned at the end of workchain execution
            """)

    def setup(self):
        """Initialize context variables that are used during the logical flow of the BaseRestartWorkChain."""
        self.ctx.unexpected_failure = False
        self.ctx.submission_failure = False
        self.ctx.restart_calc = None
        self.ctx.is_finished = False
        self.ctx.iteration = 0

        return

    def should_run_calculation(self):
        """
        Return whether a new calculation should be run.

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
        process = self._calculation_class.process()
        running = self.submit(process, **inputs)

        if hasattr(running, 'pid'):
            self.report('launching {}<{}> iteration #{}'.format(self._calculation_class.__name__, running.pid, self.ctx.iteration))
        else:
            self.report('launching {}<{}> iteration #{}'.format(self._calculation_class.__name__, running.pk, self.ctx.iteration))

        return ToContext(calculations=append_(running))

    def inspect_calculation(self):
        """
        Analyse the results of the previous calculation.

        In particular, check whether it finished successfully or if not troubleshoot the cause and handle the errors,
        or abort if unrecoverable error was found.
        """
        try:
            calculation = self.ctx.calculations[-1]
        except IndexError:
            self._fail_compat(
                exception=IndexError('the first iteration finished without returning a {}'.format(self._calculation_class.__name__)))
            return

        # Done: successful completion of last calculation
        if finished_ok_compat(calculation):
            self.report('{}<{}> completed successfully'.format(self._calculation_class.__name__, calculation.pk))
            self.ctx.restart_calc = calculation
            self.ctx.is_finished = True

        # Abort: exceeded maximum number of retries
        elif self.ctx.iteration >= self.inputs.max_iterations.value:
            self._fail_compat(
                exception=MaxIterationsFailure('reached the maximumm number of iterations {}: last ran {}<{}>'.format(
                    self.inputs.max_iterations.value, self._calculation_class.__name__, calculation.pk)))

        # Abort: unexpected state of last calculation
        elif calculation.get_state() not in self._expected_calculation_states:
            self._fail_compat(
                exception=UnexpectedCalculationFailure('unexpected state ({}) of {}<{}>'.format(
                    calculation.get_state(), self._calculation_class.__name__, calculation.pk)))

        # Retry or abort: submission failed, try to restart or abort
        elif calculation.get_state() in [calc_states.SUBMISSIONFAILED]:
            self._handle_submission_failure(calculation)
            self.ctx.submission_failure = True

        # Retry or abort: calculation finished or failed
        elif calculation.get_state() in [calc_states.FINISHED, calc_states.FAILED]:

            # Calculation was at least submitted successfully, so we reset the flag
            self.ctx.submission_failure = False

            # Check output for problems independent on calculation state and that do not trigger parser warnings
            self._handle_calculation_sanity_checks(calculation)

            # Calculation failed, try to salvage it or handle any unexpected failures
            if calculation.get_state() in [calc_states.FAILED]:
                try:
                    _ = self._handle_calculation_failure(calculation)
                except UnexpectedCalculationFailure as exception:
                    self._handle_unexpected_failure(calculation, exception)
                    self.ctx.unexpected_failure = True

            # Calculation finished: but did not finish ok, simply try to restart from this calculation
            else:
                self.ctx.unexpected_failure = False
                self.ctx.restart_calc = calculation
                self.report('calculation terminated without errors but did not complete successfully, restarting')

        return

    def results(self):
        """Attach the outputs specified in the output specification from the last completed calculation."""
        self.report('workchain completed after {} iterations'.format(self.ctx.iteration))

        for name, port in self.spec().outputs.iteritems():
            if port.required and name not in self.ctx.restart_calc.out:
                self.report('the spec specifies the output {} as required but was not an output of {}<{}>'.format(
                    name, self._calculation_class.__name__, self.ctx.restart_calc.pk))

            if name in self.ctx.restart_calc.out:
                node = self.ctx.restart_calc.out[name]
                self.out(name, self.ctx.restart_calc.out[name])
                if self._verbose:
                    self.report("attaching the node {}<{}> as '{}'".format(node.__class__.__name__, node.pk, name))

        return

    def on_destroy(self):
        """Clean remote folders of the calculations called in the workchain if the clean_workdir input is True."""
        if hasattr(super(BaseRestartWorkChain, self), 'on_destroy'):
            super(BaseRestartWorkChain, self).on_destroy()  # pylint: disable=no-member

            if not self.has_finished() or self.inputs.clean_workdir.value is False:
                return
        else:
            super(BaseRestartWorkChain, self).on_terminated()  # pylint: disable=no-member

        cleaned_calcs = []

        for calculation in self.calc.get_outputs(link_type=LinkType.CALL):
            try:
                calculation.out.remote_folder._clean()  # pylint: disable=protected-access
                cleaned_calcs.append(calculation.pk)
            except BaseException:
                pass

        if cleaned_calcs:
            self.report('cleaned remote folders of calculations: {}'.format(' '.join(map(str, cleaned_calcs))))

    # pylint: disable=no-self-use,unused-argument
    def _handle_calculation_sanity_checks(self, calculation):
        """
        Perform additional sanity checks on successfully completed calculation.

        Calculations that run successfully may still have problems that can only be determined when inspecting
        the output. The same problems may also be the hidden root of a calculation failure. For that reason,
        after verifying that the calculation ran, regardless of its calculation state, we perform some sanity
        checks.
        """
        return

    def _handle_submission_failure(self, calculation):
        """
        Handle problems during calculation submission.

        The submission of the calculation has failed. If the submission_failure flag is set to true, this
        is the second consecutive submission failure and we abort the workchain Otherwise we restart once more.
        """
        if self.ctx.submission_failure:
            self._fail_compat(
                excption=UnexpectedCalculationFailure('submission for {}<{}> failed for the second consecutive time'.format(
                    self._calculation_class.__name__, calculation.pk)))
        else:
            self.report('submission for {}<{}> failed, restarting once more'.format(self._calculation_class.__name__, calculation.pk))

        return

    def _handle_unexpected_failure(self, calculation, exception=None):
        """
        Handle unexpected failure of a calculation.

        The calculation has failed for an unknown reason and could not be handled. If the unexpected_failure
        flag is true, this is the second consecutive unexpected failure and we abort the workchain.
        Otherwise we restart once more.
        """
        if exception:
            self.report('{}'.format(exception))

        if self.ctx.unexpected_failure:
            self._fail_compat(
                exception=UnexpectedCalculationFailure('failure of {}<{}> could not be handled for the second consecutive time'.format(
                    self._calculation_class.__name__, calculation.pk)))
        else:
            self.report('failure of {}<{}> could not be handled, restarting once more'.format(self._calculation_class.__name__,
                                                                                              calculation.pk))

        return

    def _handle_calculation_failure(self, calculation):
        """
        Analyze calculation failure and prepare to restart with corrected inputs.

        The calculation has failed so we try to analyze the reason and change the inputs accordingly
        for the next calculation. If the calculation failed, but did so cleanly, we set it as the
        restart_calc, in all other cases we do not replace the restart_calc
        """
        try:
            outputs = calculation.out.output_parameters.get_dict()
            _ = outputs['warnings']
            _ = outputs['parser_warnings']
        except (AttributeError, KeyError) as exception:
            raise UnexpectedCalculationFailure(exception)

        is_handled = False

        # Sort the handlers based on their priority in reverse order
        handlers = sorted(self._error_handlers, key=lambda x: x.priority, reverse=True)

        if not handlers:
            raise UnexpectedCalculationFailure('no calculation error handlers were registered')

        for handler in handlers:
            handler_report = handler.method(self, calculation)

            # If at least one error is handled, we consider the calculation failure handled
            if handler_report and handler_report.is_handled:
                self.ctx.restart_calc = calculation
                is_handled = True

            # After certain error handlers, we may want to skip all other error handling
            if handler_report and handler_report.do_break:
                break

        # If none of the executed error handlers reported that they handled an error, the failure reason is unknown
        if not is_handled:
            raise UnexpectedCalculationFailure('calculation failure was not handled')

        return


ErrorHandler = namedtuple('ErrorHandler', 'priority method')
"""
A namedtuple to define an error handler for a :class:`~aiida.work.workchain.WorkChain`.

The priority determines in which order the error handling methods are executed, with
the higher priority being executed first. The method defines an unbound WorkChain method
that takes an instance of a :class:`~aiida.orm.implementation.general.calculation.job.AbstractJobCalculation`
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
