"""Test the Kpoints io interface"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.kpoints import KpParser


#@pytest.mark.parametrize(['vasp_kpoints'], [('list',)], indirect=True)
def test_parse_kpoints(vasp_kpoints):
    """Parse a reference KPOINTS file with the KpParser and compare the result to a reference kpoints-node."""

    kpoints, _ = vasp_kpoints

    if kpoints.get_attrs().get('mesh'):
        file_path = data_path('kpoints', 'KPOINTS_mesh')
        method = 'get_kpoints_mesh'
        param = 'mesh'
    elif kpoints.get_attrs().get('array|kpoints'):
        file_path = data_path('kpoints', 'KPOINTS_list')
        method = 'get_kpoints'
        param = 'list'

    parser = KpParser(None, path=file_path)
    result = parser.get_quantity('kpoints', {})
    if param == 'list':
        assert getattr(result['kpoints'], method)().all() == getattr(kpoints, method)().all()
    if param == 'mesh':
        assert getattr(result['kpoints'], method)() == getattr(kpoints, method)()
