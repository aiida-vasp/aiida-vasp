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
    base_parser._compose_retrieved_content()  # pylint: disable=protected-access
    assert os.path.isfile(base_parser._get_file('OUTCAR'))  # pylint: disable=protected-access
    assert os.path.exists(base_parser._get_file('OUTCAR'))  # pylint: disable=protected-access
    # This should now not eject an OSError as this is handled by the logger
    # and None is returned.
    assert base_parser._get_file('NonExistent') is None  # pylint: disable=protected-access
