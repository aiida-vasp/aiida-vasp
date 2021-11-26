"""Test the BaseFileParser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import,protected-access

import pytest

from aiida_vasp.utils.aiida_utils import get_data_node
from aiida_vasp.parsers.content_parsers.base import BaseFileParser


def test_base_object_parser(tmpdir):
    """Test functionality of the BaseFileParser, that is not already covered by any of the subclass tests."""
    with pytest.raises(TypeError):
        # Need to be initialized with either the `handler` or `data` parameters
        parser = BaseFileParser()
    temp_path = tmpdir / 'dummyfile'
    with pytest.raises(NotImplementedError):
        # Does not have a init_handler method, should be implemented in child
        with open(temp_path, 'w') as handler:
            parser = BaseFileParser(handler=handler)
    with pytest.raises(NotImplementedError):
        # Does not have a init_data method, should be implemented in child
        parser = BaseFileParser(data=get_data_node('dict', dict={}))
    with pytest.raises(TypeError):
        # Check that something else than IOBase raises an exception as it is not allowed
        parser = BaseFileParser(handler=[])
    with pytest.raises(TypeError):
        # Check that something else than Data for the AiiDA data type raises an exception.
        parser = BaseFileParser(data=[])
