"""Update version numbers everywhere based on git tags"""
import os
import re
import json
import fileinput
import contextlib
import subprocess32 as subprocess

from packaging import version


def subpath(*args):
    return os.path.realpath(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', *args))


@contextlib.contextmanager
def file_input(*args, **kwargs):
    input_fo = fileinput.FileInput(*args, **kwargs)
    try:
        yield input_fo
    finally:
        input_fo.close()


class VersionUpdater(object):
    """Version number synchronisation interface"""

    def __init__(self):
        self.top_level_init = subpath('aiida_vasp', '__init__.py')
        self.setup_json = subpath('setup.json')
        self.version = self.get_version()

    def write_to_init(self):
        with open(self.top_level_init, 'r') as init_fo:
            init_content = init_fo.read()
        with open(self.top_level_init, 'w') as init_fo:
            init_fo.write(
                re.sub(r'(__version__ = )([\'"])(.*)([\'"])',
                       r'\1\g<2>{}\4'.format(str(self.version)), init_content,
                       re.DOTALL | re.MULTILINE))

    def write_to_setup(self):
        with open(self.setup_json, 'r') as setup_fo:
            setup = json.load(setup_fo)
        setup['version'] = str(self.version)
        with open(self.setup_json, 'w') as setup_fo:
            json.dump(setup, setup_fo, indent=4, sort_keys=True)

    @staticmethod
    def get_version():
        version_byte_string = subprocess.check_output(['git', 'describe'])
        return version.parse(version_byte_string)

    def sync(self):
        self.write_to_init()
        self.write_to_setup()


if __name__ == '__main__':
    VERSION_UPDATER = VersionUpdater()
    VERSION_UPDATER.sync()
