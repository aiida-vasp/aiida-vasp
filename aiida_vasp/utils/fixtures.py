"""pytest-style test fixtures"""
# pylint: disable=unused-argument,redefined-outer-name
import pytest

from aiida.utils.fixtures import FixtureManager


@pytest.fixture(scope='session')
def aiida_env():
    manager = FixtureManager()
    manager.create_profile()
    yield manager
    manager.destroy_all()


@pytest.fixture(scope='function')
def fresh_aiida_env(aiida_env):
    yield
    aiida_env.reset_db()
