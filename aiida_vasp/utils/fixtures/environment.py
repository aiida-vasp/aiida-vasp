""" # noqa: D205
Fixtures related to environments
--------------------------------
"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
from __future__ import absolute_import
from __future__ import print_function
import pytest
from py import path as py_path  # pylint: disable=no-member,no-name-in-module

from aiida.manage.fixtures import fixture_manager


@pytest.fixture(scope='session')
def aiida_env():
    """
    Set up the db environment.

    If one has several Postgresql installations in the system, we might not pick
    up the right pg_ctl command. If this is the case, one can replace the
    call to the fixture manager below by
    fixture_manager(pgtest={'pg_ctl': '/usr/pgsql-9.6/bin/pg_ctl'})
    where it is possible to set the explicit location.
    """
    #with fixture_manager() as manager:
    with fixture_manager(pgtest={'pg_ctl': '/usr/pgsql-9.6/bin/pg_ctl'}) as manager:
        print('The root directory of the fixture manage is: {}'.format(manager.root_dir))
        # config_file = py_path.local(manager.root_dir).join('.aiida', 'config.json')
        yield manager


@pytest.fixture()
def fresh_aiida_env(aiida_env):
    """Reset the database before and after the test function."""
    yield aiida_env
    aiida_env.reset_db()
