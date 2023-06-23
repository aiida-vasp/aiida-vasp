"""Test the BaseFileParser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import,protected-access

import pytest

from aiida_vasp.parsers.content_parsers.base import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_node


def test_base_content_parser(tmpdir):
    """Test functionality of the BaseFileParser, that is not already covered by any of the subclass tests."""
    with pytest.raises(TypeError):
        # Need to be initialized with either the `handler` or `data` parameters
        _ = BaseFileParser()
    temp_path = tmpdir / 'dummyfile'
    with pytest.raises(NotImplementedError):
        # Does not have a init_handler method, should be implemented in child
        with open(temp_path, 'w', encoding='utf8') as handler:
            _ = BaseFileParser(handler=handler)
    with pytest.raises(NotImplementedError):
        # Does not have a init_data method, should be implemented in child
        _ = BaseFileParser(data=get_data_node('core.dict', dict={}))
    with pytest.raises(TypeError):
        # Check that something else than Data for the AiiDA data type raises an exception.
        _ = BaseFileParser(data=[])
