"""
Fixtures for workchains.

------------------------
Fixtures representing WorkChains, for usage in tests.
This can be either Top level workchains using mock-ups of low level WorkChains
for testing complex WorkChains, or the mock-ups themselves.
"""
# pylint: disable=too-few-public-methods, redefined-outer-name, import-outside-toplevel

import pytest


@pytest.fixture
def mock_base_workchain():
    """Fixture for a mock-up of the VaspWorkChain."""
    from aiida.plugins import WorkflowFactory
    _base_wc_cls = WorkflowFactory('vasp.vasp')

    class VaspWorkChain(_base_wc_cls):
        pass

    return VaspWorkChain


@pytest.fixture
def mock_relax_workchain(mock_base_workchain):
    """Fixture for a RelaxWorkChain using a mock-up for the lower level VaspWorkChain."""
    from aiida.plugins import WorkflowFactory
    _base_wc_cls = WorkflowFactory('vasp.relax')

    class RelaxWorkChain(_base_wc_cls):
        _base_workchain = mock_base_workchain

    return RelaxWorkChain


@pytest.fixture
def mock_verify_workchain(mock_base_workchain):  # pylint: disable=unused-argument
    """Fixture for a VerifyWorkChain using a mock-up for the lower level VaspWorkChain."""
    from aiida.plugins import WorkflowFactory
    _base_wc_cls = WorkflowFactory('vasp.verify')

    class VerifyWorkChain(_base_wc_cls):
        """An override of the verify workchain."""

        _current_iteration = 0

        @classmethod
        def define(cls, spec):
            super(_base_wc_cls, cls).define(spec)  # pylint: disable=bad-super-call
            spec.outline(cls.run)

        def run(self):  # pylint: disable=no-self-use
            # In the future, set output based on _current_iteration
            VerifyWorkChain._current_iteration += 1

    return VerifyWorkChain


@pytest.fixture
def mock_converge_workchain(mock_base_workchain):  # pylint: disable=unused-argument
    """Fixture for a ConvergenceWorkChain using a mock-up for the lower level VaspWorkChain."""
    from aiida.plugins import WorkflowFactory
    _base_wc_cls = WorkflowFactory('vasp.converge')

    class ConvergeWorkCHain(_base_wc_cls):
        _next_workchain = mock_verify_workchain

    return ConvergeWorkCHain


def mock_factory(base_class, run_method):
    """
    Return a mock-up of a class inheriting from CalcJobNode.

    :param base_class: THe base class inheriting from CalcJobNode as e.g. defined in a plugin.
    :param run_method: a method with signature run(self) that will add outputs based on the inputs to the CalcJobNode.
    """

    class MockCalcJob(base_class):
        """A mock CalcJob class."""

        run = run_method

        @classmethod
        def define(cls, spec):
            super(base_class, cls).define(spec)  # pylint: disable=bad-super-call
            spec.outline(cls.run)

    return MockCalcJob
