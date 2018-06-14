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
    parser = EigParser(file_path=path)
    bands = numpy.array([[[-1.439825, 2.964373, 2.964373, 2.964373, 7.254542, 7.254542, 7.254542, 11.451811, 11.670398, 11.670398]]])
    inputs = {}

    # occupations from bands
    vrp_bands = vasprun_parser.get_quantity('bands', {})
    inputs['occupations'] = vrp_bands['bands'].get_array('occupations')
    inputs['structure'] = vasp_structure
    result = parser.get_quantity('eigenval-bands', {}, inputs)
    assert result['eigenval-bands'].get_bands().all() == bands.all()

    # or directly
    vrp_occupations = vasprun_parser.get_quantity('occupations', {})
    inputs['occupations'] = vrp_occupations['occupations'].get_array('total')
    inputs['structure'] = vasp_structure
    result = parser.get_quantity('eigenval-bands', {}, inputs)
    assert result['eigenval-bands'].get_bands().all() == bands.all()
