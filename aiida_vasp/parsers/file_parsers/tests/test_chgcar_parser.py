"""Test the CHGCAR parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.parsers.file_parsers.chgcar import ChgcarParser


def test_parse_chgcar(aiida_env):
    """Parse a reference CHGCAR file with the ChargcarParser and compare the result to a reference string."""
    file_name = 'CHGCAR'
    path = data_path('chgcar', file_name)
    parser = ChgcarParser(file_path=path)
    result = parser.chgcar
    content = result.get_content()
    assert result.filename == file_name
    assert content == 'This is a test CHGCAR file.\n'
