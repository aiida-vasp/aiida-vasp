"""pytest-style test fixtures"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import pytest

from aiida import __version__ as aiidav
from aiida.utils.fixtures import fixture_manager


@pytest.fixture()
def aiida_version():
    """Fetches the Aiida version as a tupple."""

    # We here assume 'major.minor.bugfix' type of version.
    versions = aiidav.split('.')
    return tuple(versions)


@pytest.fixture(scope='session')
def aiida_env():
    with fixture_manager() as manager:
        yield manager


@pytest.fixture()
def fresh_aiida_env(aiida_env):
    aiida_env.reset_db()
    yield
    aiida_env.reset_db()
