"""
Unit tests for data.potcar.PotcarWalker

PotcarWalker recursively walks a directory and it's subdirectories,
searching for POTCAR files, when it encounters a tar archive, it should be extracted to
a folder on the same level as the archive and the extracted folder added to the search.
"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import pytest
from py import path as py_path  # pylint: disable=no-member,no-name-in-module

from aiida_vasp.utils.fixtures import aiida_env
from aiida_vasp.utils.fixtures.testdata import data_path


@pytest.fixture
def temp_data_folder(tmpdir):
    pot_archive = py_path.local(data_path('pot_archive.tar.gz'))
    pot_archive.copy(tmpdir)
    return tmpdir.join('pot_archive.tar.gz')


@pytest.fixture
def potcar_walker_cls(aiida_env):
    from aiida_vasp.data.potcar import PotcarWalker
    return PotcarWalker


def test_find_potcars(temp_data_folder, potcar_walker_cls):
    walker = potcar_walker_cls(temp_data_folder)
    walker.walk()
    assert len(walker.potcars) == 7
