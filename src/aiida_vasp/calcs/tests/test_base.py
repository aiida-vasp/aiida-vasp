"""Unit tests for aiida_vasp.calcs.base."""
# pylint: disable=unused-import,unused-argument,redefined-outer-name, import-outside-toplevel
import pytest
from aiida.common.extendeddicts import AttributeDict
from aiida.engine.utils import instantiate_process
from aiida.manage.manager import get_manager
from aiida_vasp.utils.fixtures import *  # pylint: disable=wildcard-import, unused-wildcard-import


def test_generate_base_calc(base_calc):
    """Test that it is possible to start the generation of an instance of the base calculation class."""

    from aiida_vasp.calcs.base import VaspCalcBase

    manager = get_manager()
    runner = manager.get_runner()
    inputs = AttributeDict()
    metadata = AttributeDict({'options': {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}})
    inputs.metadata = metadata

    # Need more input, so we will get an AttributeError because we do not pass code
    with pytest.raises(AttributeError):
        instantiate_process(runner, VaspCalcBase, **inputs)

    # Let us try again with a code attached, thus should pass
    inputs.code = base_calc.inputs.code
    instantiate_process(runner, VaspCalcBase, **inputs)
