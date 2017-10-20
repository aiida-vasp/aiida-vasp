"""pytest-style test fixtures"""
# pylint: disable=unused-argument,redefined-outer-name
import pytest

from aiida.utils.fixtures import plugin_fixture


@pytest.fixture(scope='module')
def aiida_env():
    with plugin_fixture() as manager:
        manager.create_profile()
        yield manager
        manager.destroy_all()


@pytest.fixture(scope='function')
def fresh_aiida_env(aiida_env):
    aiida_env.reset_db()
    yield
