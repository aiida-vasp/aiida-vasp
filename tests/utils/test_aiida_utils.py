"""Test aiida_utils functionss."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import,no-member

from aiida_vasp.utils.aiida_utils import (
    get_current_user,
)


def test_get_current_user(fresh_aiida_env):
    """Assert that get_current_user returns a user."""
    user = get_current_user()
    assert user.id
    assert user.first_name == ''
    assert user.last_name == ''
    assert user.email
