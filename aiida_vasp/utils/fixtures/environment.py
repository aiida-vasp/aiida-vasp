"""
Fixtures related to environments.

---------------------------------
Here we set up pytest fixtures that yield the correct AiiDA environment
to run our tests. In particular, the fresh_aiida_env is crucial in order
to set up a dummy database and configuration that can be used during
testing. AiiDA contributes with its own fixture manager, which we use.
"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
from __future__ import absolute_import
from __future__ import print_function
import pytest


@pytest.fixture()
def fresh_aiida_env(aiida_profile):
    """Reset the database before and after the test function."""
    print('The root directory of the fixture manager is: {}'.format(aiida_profile._manager.root_dir))  # pylint: disable=protected-access
    yield aiida_profile
    aiida_profile.reset_db()
