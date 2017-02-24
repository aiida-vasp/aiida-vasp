import argparse
import os
import shutil
import sys
import subprocess32 as sp

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
parser.add_argument('--patch', nargs=1, default=[None],
    metavar='AIIDA_VERSION_NR',
    help='apply patches for cli subcommands, tests and '
    'docs')
parser.add_argument(
    'aiida_root',
    help='location of the aiida root directory \
    (usually called "aiida_core" or "aiida_epfl*"')
parser.add_argument(
    'name', nargs='?', default=None, metavar='[name]',
    help='the name of the plugin package (defaults \
    to "vasp", especially helpful if multiple version \
    should be kept around.')

plugin_root = os.path.dirname(os.path.abspath(__file__))


def srcdest(paths, aiida_root, name):
    spth, dpth = paths
    src = os.path.abspath(os.path.join(plugin_root, 'aiida', spth))
    dst = os.path.abspath(
        os.path.join(aiida_root, 'aiida', os.path.dirname(dpth),
                     name or os.path.basename(dpth)))
    return src, dst


def cpdir(src, dst):
    shutil.copytree(src, dst)


def lndir(src, dst):
    os.symlink(src, dst)

install_paths = [('orm.calc.job.vasp', 'orm/calculation/job/vasp'),
                 ('parsers.plugins.vasp', 'parsers/plugins/vasp'),
                 ('orm.data.vasp', 'orm/data/vasp'),
                 ('workflows.vasp', 'workflows/vasp'),
                 ('tools.codespc.vasp', 'tools/codespecific/vasp'),
                 ('djsite.db.subtests.vasp', 'djsite/db/subtests/vasp'),
                 ('commands', 'cmdline/vasp'),
                 ('../doc/source', '../docs/source/plugins/vasp')]

colors = {'ok': '\033[92m',
          'info': '\033[91m',
          'end': '\033[0m'}

if __name__ == '__main__':
    args = parser.parse_args()
    install_cmd = args.copy and cpdir or lndir
    # install_cmd = lndir
    i = 0
    N = len(install_paths)
    for p in install_paths:
        i += 1
        src, dest = srcdest(p, os.path.abspath(args.aiida_root), args.name)
        if not os.path.exists(dest):
            print('{info}[{i}/{N}] {ok}installing:{end} {src} -> {tgt}'.format(
                src=src, tgt=dest, i=i, N=N, **colors))
        else:
            print('{info}skipping {tgt}.'.format(tgt=dest, **colors))
        if not os.path.exists(src):
            print(('{info}something went wrong with the installer: ' +
                  'source dir does not exist{end}').format(**colors))
            sys.exit()
        elif not os.path.exists(os.path.dirname(dest)):
            print(
                ('{info}it seems that the targeted aiida_root does ' +
                 'not hold an aiida distribution{end}: ' +
                 'it does not contain folder {destdir}').format(
                     destdir=os.path.dirname(dest), **colors))
            sys.exit()
        if not os.path.exists(dest):
            install_cmd(src, dest)
            print(' {ok}-> OK{end}'.format(**colors))

    if args.patch[0] == '0.5.0':
        patch_cmd = ['patch', '-d', args.aiida_root, '-b', '-p1']
        print('{ok}-> Applying CLI patch{end}'.format(**colors))
        with open(os.path.join(plugin_root, 'epfl_0.5.0.patch')) as pfile:
            try:
                sp.check_call(patch_cmd, stdin=pfile)
            except Exception as e:
                print(' {info} something went terribly wrong '
                      'while applying patches. "patch" has made backups, '
                      'but you will have to restore the originals manually. ')
                print >> stderr,  e.msg
        print(' {ok}-> OK{end}'.format(**colors))
