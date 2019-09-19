"""Unit tests for BaseParser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
import os

import pytest
from aiida.common.folders import SandboxFolder

from aiida_vasp.parsers.base import BaseParser
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import calc_with_retrieved


@pytest.fixture()
def base_parser(calc_with_retrieved):
    file_path = str(os.path.abspath(os.path.dirname(__file__)) + '/../../test_data/basic_run/')
    node = calc_with_retrieved(file_path)
    return BaseParser(node)


def test_get_file(base_parser):
    """Test getting a retrieved output file."""
    assert os.path.isfile(base_parser.get_file('OUTCAR'))
    assert os.path.exists(base_parser.get_file('OUTCAR'))
    # This should now not eject an OSError as this is handled by the logger
    # and None is returned.
    assert base_parser.get_file('NonExistent') is None
