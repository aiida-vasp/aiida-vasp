"""
Fixtures related to environments.

---------------------------------
Here we set up pytest fixtures that yield the correct AiiDA environment
to run our tests. In particular, the fresh_aiida_env is crucial in order
to set up a dummy database and configuration that can be used during
testing. AiiDA contributes with its own fixture manager, which we use.
"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name, protected-access
from __future__ import absolute_import
from __future__ import print_function
import pytest

from aiida.manage.tests import TemporaryProfileManager


@pytest.fixture()
def fresh_aiida_env(aiida_profile):
    """Reset the database before and after the test function."""
    if isinstance(aiida_profile._manager, TemporaryProfileManager):
        print('The root directory of the fixture manager is: {}'.format(aiida_profile._manager.root_dir))  # pylint: disable=protected-access
    else:
        print('Using existing profile: {}'.format(aiida_profile._manager._profile.name))
    yield aiida_profile
    aiida_profile.reset_db()
