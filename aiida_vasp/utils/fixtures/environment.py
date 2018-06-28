"""pytest-style test fixtures"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import os

import pytest
from py import path as py_path  # pylint: disable=no-member,no-name-in-module

from aiida import __version__ as aiidav
from aiida.utils.fixtures import fixture_manager


@pytest.fixture(scope='session')
def aiida_env():
    """Set up the db environment."""
    # The following should be addressed in aiida-core in the long term
    from aiida.utils import fixtures
    if hasattr(fixtures, '_GLOBAL_FIXTURE_MANAGER'):
        global_fixture_manager = fixtures._GLOBAL_FIXTURE_MANAGER  # pylint: disable=no-member,protected-access
    elif hasattr(fixtures, '_PYTEST_FIXTURE_MANAGER'):
        global_fixture_manager = fixtures._PYTEST_FIXTURE_MANAGER  # pylint: disable=no-member,protected-access
    global_fixture_manager.create_root_dir()
    root_dir = py_path.local(global_fixture_manager.root_dir)
    config_file = root_dir.join('.aiida', 'config.json')
    config_file.ensure()
    config_file.write('{}')
    os.environ['AIIDA_PATH'] = root_dir.join('.aiida').strpath
    print config_file.strpath
    with fixture_manager() as manager:
        yield manager


@pytest.fixture()
def fresh_aiida_env(aiida_env):
    aiida_env.reset_db()
    yield
    aiida_env.reset_db()
