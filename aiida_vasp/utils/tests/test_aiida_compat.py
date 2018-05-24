"""Test aiida compatibility functions."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import,no-member

import pytest

from aiida_vasp.utils.aiida_utils import get_current_user
from aiida_vasp.utils.fixtures import aiida_env


def test_get_current_user(aiida_env):
    from aiida.orm.user import User
    user = get_current_user()
    assert isinstance(user, User)
