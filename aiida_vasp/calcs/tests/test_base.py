"""Unit tests for aiida_vasp.calcs.base."""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import pytest

from aiida.common.extendeddicts import AttributeDict
from aiida.engine.utils import instantiate_process
from aiida.manage.manager import get_manager


def test_generate_base_calc(base_calc):
    """Test that it is possible to generate an instance of the base calculation class."""

    from aiida_vasp.calcs.base import VaspCalcBase

    manager = get_manager()
    runner = manager.get_runner()
    inputs = AttributeDict()

    metadata = AttributeDict({'options': {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}})

    inputs.metadata = metadata

    with pytest.raises(ValueError):
        instantiate_process(runner, VaspCalcBase, **inputs)
