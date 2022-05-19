"""
Tests for mock vasp
"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import,no-member, import-outside-toplevel, too-many-arguments
from pathlib import Path
from tempfile import mkdtemp
from shutil import rmtree, copy2

import pytest
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.mock_code import MockRegistry, MockVasp, get_hash
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_node, aiida_version, cmp_version, create_authinfo
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP


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


def test_get_hash_list():
    """Test generating has for nested dict/list"""
    dict1 = {'1': 2, 3: 4, 'a': [1, -0., '3']}
    hash1 = get_hash(dict1)[0]
    dict2 = {'1': 2, 3: 4, 'a': [1, 0., '3']}
    hash2 = get_hash(dict2)[0]

    assert hash1 == hash2

    dict1 = {'1': 2, 3: 4, 'a': [1, -0., '3', {'1': 0.}]}
    hash1 = get_hash(dict1)[0]
    dict2 = {'1': 2, 3: 4, 'a': [1, 0., '3', {'1': -0.}]}
    hash2 = get_hash(dict2)[0]

    assert hash1 == hash2


@pytest.fixture(scope='module')
def mock_registry():
    """
    Get an mock registry object
    """
    return MockRegistry()


@pytest.fixture
def custom_registry():
    """
    Return an temporary registry
    """
    temp_base = mkdtemp()
    yield MockRegistry(base_path=Path(temp_base))
    rmtree(temp_base)


@pytest.fixture
def temp_path() -> Path:
    """Return an temporary folder"""
    temp_base = mkdtemp()
    yield Path(temp_base)
    rmtree(temp_base)


def test_registry_scan(mock_registry):
    """
    Test repository scanning
    """
    mock_registry.scan()
    assert len(mock_registry.reg_hash) > 0
    # Check some existing mocks are there already
    assert 'test_bands_wc' in mock_registry.reg_name


def test_registry_extract(mock_registry):
    """Test extracting an folder from the registry"""

    tmpfolder = mkdtemp()
    mock_registry.extract_calc_by_path('test_bands_wc', tmpfolder)
    objects = [path.name for path in Path(tmpfolder).glob('*')]
    assert 'OUTCAR' in objects
    assert 'vasprun.xml' in objects
    assert 'INCAR' in objects

    rmtree(tmpfolder)


def test_registry_match(mock_registry):
    """Test round-trip hash compute and matching"""

    hash_val = mock_registry.compute_hash(mock_registry.base_path / 'test_bands_wc/inp')
    assert hash_val in mock_registry.reg_hash


def test_registry_folder_upload(mock_registry, custom_registry, temp_path):
    """Test uploading a folder to the registry"""

    # Exact an existing calculation to the folder
    mock_registry.extract_calc_by_path('test_bands_wc', temp_path)
    # Upload to a different registry
    custom_registry.upload_calc(temp_path, 'upload-example')

    # Reset the direcotry
    rmtree(str(temp_path))
    temp_path.mkdir()

    # Extract and validate
    assert 'upload-example' in custom_registry.reg_name
    custom_registry.extract_calc_by_path('upload-example', temp_path)
    objects = [path.name for path in temp_path.glob('*')]

    assert 'OUTCAR' in objects
    assert 'vasprun.xml' in objects
    assert 'INCAR' in objects


@pytest.mark.parametrize([
    'vasp_structure',
    'vasp_kpoints',
], [('str', 'mesh')], indirect=True)
def test_registry_upload_aiida(run_vasp_process, custom_registry, temp_path):
    """Test upload from an aiida calculation"""

    _, node = run_vasp_process()
    custom_registry.upload_aiida_calc(node, 'upload-example')

    # Exact the calculation
    custom_registry.extract_calc_by_path('upload-example', temp_path)

    objects = [path.name for path in temp_path.glob('*')]
    assert 'OUTCAR' in objects
    assert 'vasprun.xml' in objects
    assert 'INCAR' in objects


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_registry_upload_wc(fresh_aiida_env, run_vasp_process, custom_registry, temp_path):
    """Return an VaspWorkChain node that has been executed."""
    _, node = run_vasp_process(process_type='workchain')
    custom_registry.upload_aiida_work(node, 'upload-example')
    # Extract the calculation
    repo_path = custom_registry.get_path_by_name('upload-example/calc-000')
    objects = [path.name for path in repo_path.glob('out/*')]
    assert 'OUTCAR' in objects
    assert 'vasprun.xml' in objects

    custom_registry.extract_calc_by_path('upload-example/calc-000', temp_path)
    objects = [path.name for path in temp_path.glob('*')]
    assert 'OUTCAR' in objects
    assert 'vasprun.xml' in objects
    assert 'INCAR' in objects


def test_mock_vasp(mock_registry, temp_path):
    """Test the MockVasp class"""
    import os
    # Setup the input directory
    mock_vasp = MockVasp(temp_path, mock_registry)
    base_path = Path(data_path('test_bands_wc', 'inp'))

    with pytest.raises(ValueError):
        # Should fail due to not having any input files present and
        # we can then not match it to a test dataset
        mock_vasp.run()

    for obj in ['INCAR']:
        copy2(base_path / obj, temp_path / obj)

    with pytest.raises(ValueError):
        # Should fail due to missing POSCAR
        mock_vasp.run()

    for obj in ['POSCAR']:
        copy2(base_path / obj, temp_path / obj)

    mock_vasp.run()

    objects = [path.name for path in temp_path.glob('*')]
    assert 'OUTCAR' in objects
    assert 'vasprun.xml' in objects
