"""Test the OUTCAR io interface"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.outcar import OutcarParser


def test_parse_outcar():
    """Parse a reference OUTCAR file with the OutcarParser and compare the result to a reference value."""
    file_name = 'OUTCAR'
    path = data_path('outcar', file_name)
    parser = OutcarParser(path, file_name, None)
    params = parser.get_quantity('parameters', {})
    result = params['parameters'].get_dict()
    assert result['volume'] == 65.94
    assert result['efermi'] == 7.2948
    assert result['energies']
    assert result['symmetries']['num_space_group_operations'] == 48
    assert result['symmetries']['num_point_group_operations'] == 48
    assert result['symmetries']['point_symmetry'] == 'O_h'
    assert result['symmetries']['space_group'] == 'D_2d'
