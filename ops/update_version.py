"""Update version numbers everywhere based on git tags."""
import os
import re
import json
import fileinput
import contextlib
from pathlib import Path
import subprocess32 as subprocess

from packaging import version


def subpath(*args):
    return os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', *args))


@contextlib.contextmanager
def file_input(*args, **kwargs):
    """Context manager for a FileInput object."""
    input_fo = fileinput.FileInput(*args, **kwargs)
    try:
        yield input_fo
    finally:
        input_fo.close()


class VersionUpdater(object):  # pylint: disable=useless-object-inheritance
    """
    Version number synchronisation interface.

    Updates the version information in

    * setup.json
    * aiida_vasp/__init__.py

    to the current version number.

    The current version number is either parsed from the output of ``git describe --tags --match v*.*.*``, or if the command fails for
    any reason, from setup.json. The current version number is decided on init, syncronization can be executed by calling ``.sync()``.
    """

    version_pat = re.compile(r'\d+.\d+.\d+')
    init_version_pat = re.compile(r'(__version__ = )([\'"])(.*)([\'"])', re.DOTALL | re.MULTILINE)

    def __init__(self):
        """Initialize with documents that should be kept up to date and actual version."""
        self.top_level_init = Path(subpath('aiida_vasp', '__init__.py'))
        self.setup_json = Path(subpath('setup.json'))
        self.version = self.get_version()

    def write_to_init(self):
        init_content = self.top_level_init.read_text()
        self.top_level_init.write_text(
            re.sub(self.init_version_pat, r'\1\g<2>{}\4'.format(str(self.version)), init_content, re.DOTALL | re.MULTILINE))

    def write_to_setup(self):
        """Write the updated version number to the setup file."""
        setup = json.load(str(self.setup_json))
        setup['version'] = str(self.version)
        with open(str(self.setup_json), 'w') as setup_fo:
            json.dump(setup, setup_fo, indent=4, sort_keys=True)

    @property
    def setup_version(self):
        with self.setup_json.open() as file_object:
            return version.parse(json.load(file_object)['version'])

    @property
    def init_version(self):
        match = re.search(self.init_version_pat, self.top_level_init.read_text())
        if not match:
            raise AttributeError('No __version__ found in top-level __init__.py')
        return version.parse(match.groups()[2])

    @property
    def tag_version(self):
        """Get the current version number from ``git describe``, fall back to setup.json."""
        try:
            describe_byte_string = subprocess.check_output(  # pylint: disable=no-member
                ['git', 'describe', '--tags', '--match', 'v*.*.*'],
                universal_newlines=True,
            )
            version_string = re.findall(self.version_pat, describe_byte_string)[0]
        except subprocess.CalledProcessError:  # pylint: disable=no-member
            with open(str(self.setup_json), 'r') as setup_fo:
                setup = json.load(setup_fo)
                version_string = setup['version']

        return version.parse(version_string)

    def get_version(self):
        return max(self.setup_version, self.init_version, self.tag_version)

    def sync(self):
        if self.version > self.init_version:
            self.write_to_init()
        if self.version > self.setup_version:
            self.write_to_setup()


if __name__ == '__main__':
    VERSION_UPDATER = VersionUpdater()
    VERSION_UPDATER.sync()
