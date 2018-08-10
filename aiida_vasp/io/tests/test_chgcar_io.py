"""Test the CHGCAR io interface"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.chgcar import ChgcarParser


def test_parse_chgcar():
    """Parse a reference CHGCAR file with the ChargcarParser and compare the result to a reference string."""
    file_name = 'CHGCAR'
    path = data_path('chgcar', file_name)
    parser = ChgcarParser(file_path=path)
    result = parser.chgcar
    with open(result.get_file_abs_path(), 'r') as file_obj:
        content = file_obj.readline()
    assert result.filename == file_name
    assert content == 'This is a test CHGCAR file.\n'
