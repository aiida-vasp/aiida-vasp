"""Unit tests for aiida_vasp.calcs.base"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import pytest
from aiida.common import ValidationError
from aiida.engine.utils import instantiate_process
from aiida.orm import Code
from aiida.manage.manager import get_manager
from aiida_vasp.utils.fixtures import *


@pytest.fixture()
def base_calc(fresh_aiida_env):
    from aiida_vasp.calcs.base import VaspCalcBase
    manager = get_manager()
    runner = manager.get_runner()
    instantiate_process(runner, VaspCalcBase, code=Code())


def test_generate_base_calc(base_calc):
    """Test that it is possible to generate an instance of the base calculation class."""
    new_calc = base_calc()
    print(new_calc)
    assert False
