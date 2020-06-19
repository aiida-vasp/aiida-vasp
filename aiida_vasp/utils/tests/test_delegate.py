# pylint: disable=unused-import,too-few-public-methods,missing-docstring,no-self-use,no-member
"""Test the Delegate class."""

import pytest

from aiida_vasp.utils.delegates import Delegate


class VaspParser(object):  # pylint: disable=useless-object-inheritance

    def __init__(self):
        super(VaspParser, self).__init__()

        # Initialise the 'get_quantity' delegate:
        setattr(self, 'get_quantity', Delegate())


class FileParser(object):  # pylint: disable=useless-object-inheritance

    def __init__(self, cls):
        self._vasp_parser = cls
        if cls is not None:
            self._vasp_parser.get_quantity.append(self.parse_file)

    def parse_file(self, inputs):
        return inputs


def test_delegate():
    """Test the functionality of the Delegate class wrapper."""

    vasp_parser = VaspParser()
    assert not vasp_parser.get_quantity
    file_parser = FileParser(vasp_parser)
    inputs = 'This is a test input'
    assert vasp_parser.get_quantity(inputs) == inputs
    vasp_parser.get_quantity.remove(file_parser.parse_file)
    assert vasp_parser.get_quantity == []
