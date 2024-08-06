"""Test the get_bands_info method"""

# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import,no-member, import-outside-toplevel
import numpy as np
import pytest
from aiida_vasp.utils.compare_bands import (
    get_band_properties,
    get_band_properties_from_data,
)


@pytest.fixture
def example_bands(fresh_aiida_env):
    """Example eigen values and occupations"""
    from aiida.orm import BandsData

    bdata = BandsData()
    bdata.set_kpoints([[0.0, 0.0, 0.0]])
    occ = np.array([[1, 1, 1, 0, 0, 0]], dtype=np.float64)
    eigen = np.array([[-1, -1, -1, 0, 0, 0]])
    bdata.set_bands(bands=eigen, occupations=occ)
    return bdata


@pytest.fixture
def example_bands_v2(fresh_aiida_env):
    """Example eigen values and occupations"""
    from aiida.orm import BandsData

    bdata = BandsData()
    bdata.set_kpoints([[0, 0, 0], [0.25, 0.25, 0.25]])
    occ = np.array([[1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0]], dtype=np.float64)
    eigen = np.array([[-1, -1, -1, 0, 0, 0], [-0.5, -0.5, -0.5, 0, 0, 0]])
    bdata.set_bands(bands=eigen, occupations=occ)
    return bdata


def test_bands_info():
    """Test getting bands info"""
    occ = np.array([[1, 1, 1, 0, 0, 0]])
    eigen = np.array([[-1, -1, -1, 0, 0, 0]])
    bands_info = get_band_properties(eigen, occ)
    assert bands_info['is_direct_gap'] is True
    assert bands_info['band_gap'] == pytest.approx(1.0)

    occ = np.array([[1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0]])
    eigen = np.array([[-1, -1, -1, 0, 0, 0], [-0.5, -0.5, -0.5, 0, 0, 0]])
    bands_info = get_band_properties(eigen, occ)
    assert bands_info['is_direct_gap'] is False
    assert bands_info['band_gap'] == pytest.approx(0.5)


def test_bands_info_from_data(example_bands, example_bands_v2):
    """Test getting info from BandsData"""
    bands_info = get_band_properties_from_data(example_bands)
    assert bands_info['is_direct_gap'] is True
    assert bands_info['band_gap'] == pytest.approx(1.0)

    bands_info = get_band_properties_from_data(example_bands_v2)
    assert bands_info['is_direct_gap'] is False
    assert bands_info['band_gap'] == pytest.approx(0.5)
    assert bands_info['vbm'] == pytest.approx(-0.5)
    assert bands_info['cbm'] == pytest.approx(0.0)
