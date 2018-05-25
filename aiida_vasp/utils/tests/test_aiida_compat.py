"""Test aiida compatibility functions."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import,no-member

import pytest

from aiida_vasp.utils.aiida_utils import get_current_user, backend_obj_users
from aiida_vasp.utils.fixtures import aiida_env


def test_get_current_user(aiida_env):
    """Assert that get_current_user returns a user in all tested aiida versions."""
    from aiida.orm.user import User
    user = get_current_user()
    if backend_obj_users():
        assert isinstance(user, User)
    assert user.first_name
    assert user.last_name
    assert user.email
