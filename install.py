import argparse
import os
import shutil
import sys

parser = argparse.ArgumentParser(
    description='install the plugin files into the \
    given aiida folder, defaults to symlinking but \
    copy can be used instead.')
parser.add_argument(
    '--copy',
    action='store_true',
    help='instead of symlinking into aiida, copy \
    the plugin files. helpful during development \
    to keep a stable version around')
parser.add_argument(
    'aiida_root',
    help='location of the aiida root directory \
    (usually called "aiida_core" or "aiida_epfl*"')
parser.add_argument(
    'name', default=None, metavar='[name]',
    help='the name of the plugin package (defaults \
    to "vasp", especially helpful if multiple version \
    should be kept around.')

plugin_root = os.path.dirname(os.path.abspath(__file__))


def srcdest(path, aiida_root, name):
    src = os.path.join(plugin_root, 'aiida', path)
    dst = os.path.join(aiida_root, 'aiida', os.path.dirname(path),
                       name or os.path.basename(path))
    return src, dst


def cpdir(src, dst):
    shutil.copytree(src, dst)


def lndir(src, dst):
    os.symlink(src, dst)

install_paths = ['orm/calculation/job/vasp',
                 'parsers/plugins/vasp',
                 'workflows/vasp']

colors = {'ok': '\033[92m',
          'info': '\033[91m',
          'end': '\033[0m'}

if __name__ == '__main__':
    args = parser.parse_args()
    install_cmd = args.copy and cpdir or lndir
    i = 0
    N = len(install_paths)
    for p in install_paths:
        i += 1
        src, dest = srcdest(p, os.path.abspath(args.aiida_root), args.name)
        print('{info}[{i}/{N}] {ok}installing:{end} {src} -> {tgt}'.format(
            src=src, tgt=dest, i=i, N=N, **colors))
        if not os.path.exists(src):
            print('{info}something went wrong with the installer: ' +
                  'source dir does not exist{end}'.format(**colors))
            sys.exit()
        elif not os.path.exists(os.path.dirname(dest)):
            print('{info}it seems that the targeted aiida_root does ' +
                  'not hold an aiida distribution{end}'.format(**colors))
            sys.exit()
        install_cmd(src, dest)
        print(' {ok}-> OK{end}'.format(**colors))
