"""Test the WAVECAR parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.parsers.file_parsers.wavecar import WavecarParser


def test_parse_wavecar(fresh_aiida_env):
    """Parse a reference CHGCAR file with the ChargcarParser and compare the result to a reference string."""
    path = data_path('wavecar', 'WAVECAR')
    parser = WavecarParser(file_path=path)
    result = parser.wavecar
    content = result.get_content()
    assert result.filename == 'WAVECAR'
    assert content == 'This is a test WAVECAR file.\n'
