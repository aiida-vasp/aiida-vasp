"""Test the EIGENVAL io interface"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.eigenval import EigParser


@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_parse_eigenval(vasprun_parser, vasp_structure):
    """Parse a reference EIGENVAL file with the EigParser and compare the result to a reference."""
    file_name = 'EIGENVAL'
    path = data_path('eigenval', file_name)
    parser = EigParser(path, file_name, None)
    bands = numpy.array([[[-1.439825, 2.964373, 2.964373, 2.964373, 7.254542, 7.254542, 7.254542, 11.451811, 11.670398, 11.670398]]])
    inputs = {}
    inputs['occupations'] = vasprun_parser.occupations
    inputs['structure'] = vasp_structure
    result = parser.get_quantity(None, 'bands', {}, inputs)
    assert result['bands'].get_bands().all() == bands.all()
