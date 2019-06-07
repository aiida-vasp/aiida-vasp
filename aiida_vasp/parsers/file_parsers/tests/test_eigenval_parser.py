"""Test the EIGENVAL parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.parsers.file_parsers.eigenval import EigParser


def test_parse_eigenval():
    """Parse a reference EIGENVAL file with the EigParser and compare the result to a reference."""
    file_name = 'EIGENVAL'
    path = data_path('eigenval', file_name)
    parser = EigParser(file_path=path)
    bands = numpy.array([[[-1.439825, 2.964373, 2.964373, 2.964373, 7.254542, 7.254542, 7.254542, 11.451811, 11.670398, 11.670398]]])
    inputs = {}

    result = parser.get_quantity('eigenval-eigenvalues', inputs)
    assert result['eigenval-eigenvalues'].all() == bands.all()
