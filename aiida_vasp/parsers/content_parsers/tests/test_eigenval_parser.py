"""Test the EIGENVAL parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.parsers.content_parsers.eigenval import EigParser


def test_parse_eigenval():
    """Parse a reference EIGENVAL with the EigParser and compare the result to a reference."""
    name = 'EIGENVAL'
    path = data_path('eigenval', name)
    parser = EigParser(path=path)
    bands = numpy.array([[[-1.439825, 2.964373, 2.964373, 2.964373, 7.254542, 7.254542, 7.254542, 11.451811, 11.670398, 11.670398]]])
    inputs = {}

    result = parser.get_quantity_from_inputs('eigenval-eigenvalues', inputs, None)
    assert result.all() == bands.all()
