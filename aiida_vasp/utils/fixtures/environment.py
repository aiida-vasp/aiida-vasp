"""pytest-style test fixtures"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import pytest

from aiida.utils.fixtures import fixture_manager

@pytest.fixture()
def element_x_not_present():
    # Element X was added in the same version as
    # the get_strict_version function.
    try:
        from aiida import get_strict_version
        return False
    except ImportError:
        return True

@pytest.fixture(scope='session')
def aiida_env():
    with fixture_manager() as manager:
        yield manager


@pytest.fixture()
def fresh_aiida_env(aiida_env):
    aiida_env.reset_db()
    yield
    aiida_env.reset_db()
