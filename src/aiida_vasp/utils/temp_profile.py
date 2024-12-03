"""
Code for getting a temporary profile for testing purposes.
"""

from aiida import load_profile, manage, orm
from aiida.storage.sqlite_temp import SqliteTempBackend

from aiida_vasp.utils.mock_code import VaspMockRegistry

__all__ = ('VaspMockRegistry', 'load_temp_profile', 'orm')


def load_temp_profile():
    """Load a temporary profile for testing/demo purposes."""
    profile = load_profile(
        SqliteTempBackend.create_profile('myprofile', options={'runner.poll.interval': 1}, debug=False),
        allow_switch=True,
    )
    config = manage.get_config()
    config.add_profile(profile)
    # Enable caching
    config.set_option('caching.enabled_for', ['aiida.calculations:vasp.vasp'])
    return profile
