"""
Fixtures representing WorkChains, for usage in tests.

This can be either Top level workchains using mock-ups of low level WorkChains
for testing complex WorkChains, or the mock-ups themselves.
"""
# pylint: disable=too-few-public-methods, redefined-outer-name

import pytest


@pytest.fixture
def mock_base_wc():
    """Fixture for a mock-up of the VaspBaseWf."""
    from aiida.orm import WorkflowFactory
    _base_wc_cls = WorkflowFactory('vasp.vasp')

    class VaspWorkChain(_base_wc_cls):
        pass

    return VaspWorkChain


@pytest.fixture
def mock_relax_wc(mock_base_wc):
    """Fixture for a VaspRelaxWf using a mock-up for the lower level BaseWf."""
    from aiida.orm import WorkflowFactory
    _base_wc_cls = WorkflowFactory('vasp.relax')

    class RelaxWorkChain(_base_wc_cls):
        _base_workchain = mock_base_wc

    return RelaxWorkChain
