"""
Code for getting a temporary profile for testing purposes.
"""

import os
from pathlib import Path
from subprocess import check_output

from aiida import load_profile, manage, orm
from aiida.manage import Profile
from aiida.storage.sqlite_temp import SqliteTempBackend

from aiida_vasp.utils.mock_code import VaspMockRegistry

__all__ = ('VaspMockRegistry', 'load_temp_profile', 'orm')


def load_temp_profile(force: bool = True) -> Profile:
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


def load_temp_profile_with_mock(force: bool = True) -> Profile:
    """Load a temporary profile with mock VASP codes for testing/demo purposes."""
    # Skip creation
    profile = load_temp_profile()

    # Register mock codes
    comp = orm.Computer('localhost', 'localhost', transport_type='core.local', scheduler_type='core.direct')
    comp.store()
    comp.set_workdir('/tmp/aiida_run/')
    comp.configure()

    vasp_path = check_output(['which', 'mock-vasp'], universal_newlines=True).strip()
    vasp_code = orm.InstalledCode(comp, vasp_path[0], default_calc_job_plugin='vasp.vasp')
    print(vasp_path[0])
    vasp_code.label = 'vasp'
    vasp_code.store()

    vasp_code = orm.InstalledCode(comp, vasp_path[0], default_calc_job_plugin='vasp.vasp')
    print(vasp_path[0])
    vasp_code.label = 'mock-vasp'
    vasp_code.store()

    os.environ['MOCK_VASP_REG_BASE'] = str((Path() / 'mock_registry').absolute())
    os.environ['MOCK_VASP_UPLOAD_PREFIX'] = 'singlepoint'

    return profile
