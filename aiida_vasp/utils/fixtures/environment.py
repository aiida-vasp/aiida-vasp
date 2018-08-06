"""pytest-style test fixtures"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import os
import json

import pytest
from py import path as py_path  # pylint: disable=no-member,no-name-in-module

from aiida import __version__ as aiidav
from aiida.utils.fixtures import fixture_manager


@pytest.fixture(scope='session')
def aiida_env():
    """Set up the db environment."""
    with fixture_manager() as manager:
        config_file = py_path.local(manager.root_dir).join('.aiida', 'config.json')
        print config_file.read()
        yield manager


@pytest.fixture()
def fresh_aiida_env(aiida_env):
    """Reset the database before and after the test function."""
    aiida_env.reset_db()
    yield aiida_env
    aiida_env.reset_db()
