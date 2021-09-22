"""Test the CHGCAR parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.parsers.object_parsers.chgcar import ChgcarParser


def test_parse_chgcar(fresh_aiida_env):
    """Parse a reference CHGCAR with the ChargcarParser and compare the result to a reference string."""
    name = 'CHGCAR'
    path = data_path('chgcar', name)
    parser = ChgcarParser(path=path)
    result = parser.chgcar
    content = result.get_content()
    assert result.filename == name
    assert content == 'This is a test CHGCAR.\n'
