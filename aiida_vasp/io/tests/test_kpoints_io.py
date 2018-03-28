"""Test the Kpoints io interface"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.fixtures.data import vasp_kpoints
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.kpoints import KpParser


def test_parse_kpoints(vasp_kpoints):
    """Parse a reference KPOINTS file with the KpParser and compare the result to a reference kpoints-node."""
    kpoints_path = data_path('kpoints', 'KPOINTS')
    kpp = KpParser(kpoints_path, 'KPOINTS')
    result = {}
    kpp.get_quantity('kpoints_raw', result)
    _, reference = vasp_kpoints
    assert result['kpoints_raw'] == reference
