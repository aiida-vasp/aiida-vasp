"""Test the Kpoints io interface"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.wavecar import WavecarParser


def test_parse_wavecar():
    """Parse a reference CHGCAR file with the ChargcarParser and compare the result to a reference string."""
    path = data_path('wavecar', 'WAVECAR')
    parser = WavecarParser(path, 'WAVECAR')
    result = {}
    parser.get_quantity('wavecar', result)
    with open(result['wavecar'].get_file_abs_path(), 'r') as file_obj:
        content = file_obj.readline()
    assert result['wavecar'].filename == 'WAVECAR'
    assert content == 'This is a test WAVECAR file.\n'
