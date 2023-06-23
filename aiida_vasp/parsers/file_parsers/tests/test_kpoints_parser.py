"""Test the KPOINTS parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.parsers.file_parsers.kpoints import KpointsParser


#@pytest.mark.parametrize(['vasp_kpoints'], [('list',)], indirect=True)
def test_parse_kpoints(vasp_kpoints):
    """
    Parse a reference KPOINTS file.

    Using the KpointsParser and compare the result to a reference
    kpoints-node.

    """

    kpoints, _ = vasp_kpoints

    try:
        _ = kpoints.get_attribute('mesh')
        file_path = data_path('kpoints', 'KPOINTS_mesh')
        method = 'get_kpoints_mesh'
        param = 'mesh'
    except AttributeError:
        pass

    try:
        _ = kpoints.get_attribute('array|kpoints')
        file_path = data_path('kpoints', 'KPOINTS_list')
        method = 'get_kpoints'
        param = 'list'
    except AttributeError:
        pass

    parser = KpointsParser(file_path=file_path)
    result = parser.kpoints
    if param == 'list':
        assert getattr(result, method)().all() == getattr(kpoints, method)().all()
    if param == 'mesh':
        assert getattr(result, method)() == getattr(kpoints, method)()


def test_parse_kpoints_write(vasp_kpoints, tmpdir):
    """
    Parse a reference KPOINTS file.

    Using the KpointsParser and compare the result to a reference
    kpoints-node.

    """

    kpoints, _ = vasp_kpoints

    try:
        _ = kpoints.get_attribute('mesh')
        method = 'get_kpoints_mesh'
        param = 'mesh'
    except AttributeError:
        pass

    try:
        _ = kpoints.get_attribute('array|kpoints')
        method = 'get_kpoints'
        param = 'list'
    except AttributeError:
        pass

    parser = KpointsParser(data=kpoints)
    temp_file = str(tmpdir.join('KPOINTS'))
    parser.write(file_path=temp_file)
    parser_reparse = KpointsParser(file_path=temp_file)
    result = parser_reparse.kpoints
    if param == 'list':
        assert getattr(result, method)().all() == getattr(kpoints, method)().all()
    if param == 'mesh':
        assert getattr(result, method)() == getattr(kpoints, method)()
