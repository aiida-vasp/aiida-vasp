"""Test the BaseFileParser"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import,protected-access

import pytest

from aiida_vasp.io.parser import BaseFileParser


def test_base_file_parser():
    """Test functionality of the BaseFileParser, that is not already covered by any of the subclass tests."""
    parser = BaseFileParser(None)
    inputs = {}
    assert parser.parsable_items == {}
    assert parser._parsed_data == {}
    with pytest.raises(Exception):
        parser._parse_file(inputs)
