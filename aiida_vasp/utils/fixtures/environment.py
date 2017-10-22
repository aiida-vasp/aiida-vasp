"""pytest-style test fixtures"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import pytest

from aiida.utils.fixtures import fixture_manager


@pytest.fixture(scope='session')
def aiida_env():
    with fixture_manager() as manager:
        yield manager


@pytest.fixture()
def fresh_aiida_env(aiida_env):
    yield
    aiida_env.reset_db()
