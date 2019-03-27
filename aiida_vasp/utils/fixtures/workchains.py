"""
Fixtures representing WorkChains, for usage in tests.

This can be either Top level workchains using mock-ups of low level WorkChains
for testing complex WorkChains, or the mock-ups themselves.
"""
# pylint: disable=too-few-public-methods, redefined-outer-name

import pytest


@pytest.fixture
def mock_base_wc():
    """Fixture for a mock-up of the VaspWorkChain."""
    from aiida.orm import WorkflowFactory
    _base_wc_cls = WorkflowFactory('vasp.vasp')

    class VaspWorkChain(_base_wc_cls):
        pass

    return VaspWorkChain


@pytest.fixture
def mock_relax_wc(mock_base_wc):
    """Fixture for a RelaxWorkChain using a mock-up for the lower level VaspWorkChain."""
    from aiida.orm import WorkflowFactory
    _base_wc_cls = WorkflowFactory('vasp.relax')

    class RelaxWorkChain(_base_wc_cls):
        _base_workchain = mock_base_wc

    return RelaxWorkChain


@pytest.fixture
def mock_verify_workchain():
    """Fixture for a VerifyWorkChain using a mock-up for the lower level VaspWorkChain."""
    from aiida.orm import WorkflowFactory
    _base_wc_cls = WorkflowFactory('vasp.verify')

    class VerifyWorkChain(_base_wc_cls):

        _current_iteration = 0

        @classmethod
        def define(cls, spec):
            super(_base_wc_cls, cls).define(spec)
            spec.outline(cls.run)

        def run(self):
            # TODO Set output based on _current_iteration
            VerifyWorkChain._current_iteration += 1

    return VerifyWorkChain


@pytest.fixture
def mock_converge_workchain(mock_verify_workchain):
    from aiida.orm import WorkflowFactory
    _base_wc_cls = WorkflowFactory('vasp.converge')

    class ConvergeWorkCHain(_base_wc_cls):
        _next_workchain = mock_verify_workchain

    return ConvergeWorkCHain
