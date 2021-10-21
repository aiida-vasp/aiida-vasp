"""Unit tests for BaseParser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
import os
import io

import pytest
from aiida.common.folders import SandboxFolder

from aiida_vasp.parsers.base import BaseParser
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import calc_with_retrieved


@pytest.fixture()
def base_parser(calc_with_retrieved):
    path = str(os.path.abspath(os.path.dirname(__file__)) + '/../../test_data/basic_run/')
    node = calc_with_retrieved(path)
    return BaseParser(node)


def test_get_object(base_parser):
    """Test getting a retrieved output object."""
    base_parser._compose_retrieved_content()  # pylint: disable=protected-access
    with base_parser._get_handler('OUTCAR') as handler:  # pylint: disable=protected-access
        # Check that it is a the right type of object.
        assert isinstance(handler, io.TextIOWrapper)
        # Check that we can do file like operations on the handlers
        content = handler.readlines()
        assert content
        # Do basic check that we can indeed access the known content
        assert content[0] == ' vasp.5.3.5 31Mar14 (build Jul 23 2014 14:06:12) complex\n'
    with base_parser._get_handler('NotExistent') as handler:  # pylint: disable=protected-access
        assert handler is None
