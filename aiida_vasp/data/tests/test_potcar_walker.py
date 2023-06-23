"""
Unit tests for data.potcar.PotcarWalker.

PotcarWalker recursively walks a directory and it's subdirectories,
searching for POTCAR files, when it encounters a tar archive, it should be extracted to
a folder on the same level as the archive and the extracted folder added to the search.
"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name, import-outside-toplevel
import shutil
from pathlib import Path
import pytest

from aiida_vasp.utils.fixtures import fresh_aiida_env
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.fixtures.data import temp_pot_folder


@pytest.fixture
def temp_data_folder(tmpdir):
    pot_archive = Path(data_path('pot_archive.tar.gz'))
    shutil.copy(pot_archive, tmpdir)
    return tmpdir.join('pot_archive.tar.gz')


@pytest.fixture
def potcar_walker_cls(fresh_aiida_env):
    from aiida_vasp.data.potcar import PotcarWalker
    return PotcarWalker


def test_find_potcars(potcar_walker_cls, temp_data_folder):
    """Make sure the walker finds the right number fo POTCAR files."""
    potcar_archive = Path(data_path('.')) / 'pot_archive'
    walker = potcar_walker_cls(str(temp_data_folder))
    walker.walk()
    assert len(walker.potcars) == 7
    assert not potcar_archive.exists()
