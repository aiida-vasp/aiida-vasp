"""Test the DOSCAR io interface"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.doscar import DosParser


def test_parse_doscar(vasprun_parser):
    """Parse a reference DOSCAR file with the DosParser and compare the result to a reference."""
    file_name = 'DOSCAR'
    path = data_path('doscar', file_name)
    parser = DosParser(path, file_name, None)
    dos = numpy.array(
        [(-3.4398, -1.10400000e-43, -2.09900000e-43), (-1.5387, 1.40000000e-01, 2.66100000e-01), (0.3624, -3.62400000e-73, 2.00000000e+00),
         (2.2636, -1.33800000e-05, 2.00000000e+00), (4.1647, 3.15600000e+00, 8.00000000e+00), (6.0659, -2.41200000e-15, 8.00000000e+00),
         (7.967, 3.15600000e+00, 1.40000000e+01), (9.8681, -1.38100000e-27, 1.40000000e+01), (11.7693, 2.90100000e+00, 1.95152000e+01),
         (13.6704, 0.00000000e+00, 2.00000000e+01)],
        dtype=[('energy', '<f8'), ('total', '<f8'), ('integrated', '<f8')])
    inputs = {}
    inputs['vrp_pdos'] = vasprun_parser.vrp_pdos
    inputs['vrp_tdos'] = vasprun_parser.vrp_tdos
    result = parser.get_quantity(None, 'dos', {}, inputs)
    # The 'tdos' array is a nested array for some reason.
    result_dos = result['dos'].get_array('tdos')[0]
    for i in range(0, dos.size):
        assert result_dos[i] == dos[i]
