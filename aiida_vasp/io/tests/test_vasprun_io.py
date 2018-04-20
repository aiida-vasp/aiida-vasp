"""Test the vasprun.xml io interface"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.vasprun import VasprunParser


def test_parse_vasprun():
    """Parse a reference vasprun.xml file with the VasprunParser and compare the result to a reference string."""
    file_name = 'vasprun.xml'
    path = data_path('vasprun', file_name)
    parser = VasprunParser(path=path)
    occupations = numpy.array([[[1., 1., 1., 1., 0.6667, 0.6667, 0.6667, -0., -0., -0.]]])
    result = parser.get_quantity('occupations', {})
    assert result['occupations'].all() == occupations.all()
    assert parser.efermi == 7.29482275
