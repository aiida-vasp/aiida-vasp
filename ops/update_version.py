"""Update version numbers everywhere based on git tags."""
import os
import re
import json
import fileinput
import contextlib
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


class VersionUpdater(object):
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

    def __init__(self):
        """Initialize with documents that should be kept up to date and actual version."""
        self.top_level_init = subpath('aiida_vasp', '__init__.py')
        self.setup_json = subpath('setup.json')
        self.version = self.get_version()

    def write_to_init(self):
        with open(self.top_level_init, 'r') as init_fo:
            init_content = init_fo.read()
        with open(self.top_level_init, 'w') as init_fo:
            init_fo.write(
                re.sub(r'(__version__ = )([\'"])(.*)([\'"])', r'\1\g<2>{}\4'.format(str(self.version)), init_content,
                       re.DOTALL | re.MULTILINE))

    def write_to_setup(self):
        """Write the updated version number to the setup file."""
        with open(self.setup_json, 'r') as setup_fo:
            setup = json.load(setup_fo)
        setup['version'] = str(self.version)
        with open(self.setup_json, 'w') as setup_fo:
            json.dump(setup, setup_fo, indent=4, sort_keys=True)

    def get_version(self):
        """Get the current version number from ``git describe``, fall back to setup.json."""
        try:
            describe_byte_string = subprocess.check_output(['git', 'describe', '--tags', '--match', 'v*.*.*'])
            version_string = re.findall(self.version_pat, describe_byte_string)[0]
        except subprocess.CalledProcessError:
            with open(self.setup_json, 'r') as setup_fo:
                setup = json.load(setup_fo)
                version_string = setup['version']

        return version.parse(version_string)

    def sync(self):
        self.write_to_init()
        self.write_to_setup()


if __name__ == '__main__':
    VERSION_UPDATER = VersionUpdater()
    VERSION_UPDATER.sync()
