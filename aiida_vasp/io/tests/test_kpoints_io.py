"""Test the Kpoints io interface"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.kpoints import KpParser


@pytest.mark.parametrize(['vasp_kpoints'], [('list',)], indirect=True)
def test_parse_kpoints(vasp_kpoints):
    """Parse a reference KPOINTS file with the KpParser and compare the result to a reference kpoints-node."""
    file_path = data_path('kpoints', 'KPOINTS')
    parser = KpParser(file_path)
    result = parser.get_quantity('kpoints', {})
    kpoints, _ = vasp_kpoints
    assert result['kpoints'].get_kpoints().all() == kpoints.get_kpoints().all()
