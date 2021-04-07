"""
Tests for mock vasp
"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import,no-member, import-outside-toplevel
from pathlib import Path
from tempfile import mkdtemp
from shutil import rmtree

import pytest
from aiida_vasp.utils.mock_vasp import MockRegistry, MockVasp, get_hash


def test_get_hash():
    """
    Test the get_hash function
    """
    dict1 = {'1': 2, 3: 4, 'a': [1, '2', '3']}
    hash1 = get_hash(dict1)[0]
    dict2 = {'1': 38, 3: 4, 'a': [1, '2', '3']}
    hash2 = get_hash(dict2)[0]
    dict3 = {3: 4, '1': 2, 'a': [1, '2', '3']}
    hash3 = get_hash(dict3)[0]

    assert hash1 != hash2
    assert hash1 == hash3


@pytest.fixture(scope='module')
def mock_registry():
    """
    Get an mock registry object
    """
    return MockRegistry()


def test_registry_scan(mock_registry):
    """
    Test repository scaning
    """
    mock_registry.scan()
    assert len(mock_registry.reg_hash) > 0
    # Check some existing mocks are there already
    assert 'test_bands_wc' in mock_registry.reg_name


def test_registry_extract(mock_registry):
    """Test extracting an folder from the registry"""

    tmpfolder = mkdtemp()
    mock_registry.extract_calc_by_path('test_bands_wc', tmpfolder)
    files = [path.name for path in Path(tmpfolder).glob('*')]
    assert 'OUTCAR' in files
    assert 'vasprun.xml' in files
    assert 'INCAR' in files

    rmtree(tmpfolder)
