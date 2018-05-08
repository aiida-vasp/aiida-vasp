"""Test the delegate decorator"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import,no-member,no-self-use,missing-docstring

import pytest

from aiida_vasp.utils.delegates import delegate


class VaspParser(object):

    @delegate()
    def get_quantity(self, inputs):
        return None


class FileParser(object):

    def __init__(self, cls):
        self._vasp_parser = cls
        if cls is not None:
            self._vasp_parser.get_quantity.add_listener(self.parse_file)

    def parse_file(self, inputs):
        return inputs


def test_delegate():
    """Test the functionality of the delegate() wrapper"""

    vasp_parser = VaspParser()
    assert not vasp_parser.get_quantity.listeners
    file_parser = FileParser(vasp_parser)
    inputs = 'This is a test input'
    assert vasp_parser.get_quantity(inputs) == inputs
    vasp_parser.get_quantity.remove_listener(file_parser.parse_file)
    assert vasp_parser.get_quantity.listeners == []
