"""
Fixtures representing WorkChains, for usage in tests.

This can be either Top level workchains using mock-ups of low level WorkChains
for testing complex WorkChains, or the mock-ups themselves.
"""
# pylint: disable=too-few-public-methods, redefined-outer-name

import pytest


@pytest.fixture
def mock_base_wf():
    """Fixture for a mock-up of the VaspBaseWf."""
    from aiida.orm import WorkflowFactory
    _base_wf_cls = WorkflowFactory('vasp.base')

    class BaseWorkChain(_base_wf_cls):
        pass

    return BaseWorkChain


@pytest.fixture
def mock_relax_wf(mock_base_wf):
    """Fixture for a VaspRelaxWf using a mock-up for the lower level BaseWf."""
    from aiida.orm import WorkflowFactory
    _base_wf_cls = WorkflowFactory('vasp.relax')

    class RelaxWorkChain(_base_wf_cls):
        _base_workchain = mock_base_wf

    return RelaxWorkChain
