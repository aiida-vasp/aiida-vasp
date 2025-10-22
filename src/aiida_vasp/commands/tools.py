"""
Module for easy access to aiida-vasp calculations.
"""

import os
import subprocess
import sys
import tempfile
from functools import wraps
from shutil import copyfileobj
from typing import Any, Callable

import click
from aiida import orm
from aiida.cmdline.commands.cmd_calcjob import calcjob_outputcat
from aiida.cmdline.params.arguments import (
    CALCULATION,
    PROCESS,
    WORKFLOW,
)
from aiida.cmdline.utils.echo import (
    echo_critical,
    echo_error,
    echo_info,
    echo_success,
)
from aiida.common.escaping import escape_for_bash
from aiida.common.exceptions import NotExistent
from aiida.orm import CalcJobNode, QueryBuilder, WorkChainNode

from aiida_vasp.utils.aiida_utils import cmp_load_verdi_data
from aiida_vasp.utils.export import export_vasp

VERDI_DATA = cmp_load_verdi_data()


@VERDI_DATA.group('vasp-tools')
def tools() -> None:
    """Tool for aiida-vasp related data"""


@tools.command('export')
@PROCESS('process')
@click.argument('folder')
@click.option(
    '--include-potcar',
    default=False,
    is_flag=True,
    help='Whether to include POTCAR in the export folder',
)
@click.option(
    '--decompress',
    default=False,
    is_flag=True,
    help='Wether to decompress the contents',
)
def export(process: orm.ProcessNode, folder: str, decompress: bool, include_potcar: bool) -> None:
    """Export a VASP calculation, works for both `VaspCalculation` or `VaspWorkChain`"""
    export_vasp(process, folder, decompress, include_potcar)


def select_calcjob_from_work(func: Callable) -> Callable:
    """Select calcjob from work"""

    @wraps(func)
    def wrapper(*args, **kwargs):
        calcjob = kwargs['calcjob']
        index = kwargs.pop('index', 0)

        if isinstance(calcjob, orm.WorkChainNode):
            valid = [process for process in calcjob.called_descendants if not process.is_finished]
            if len(valid) > 0:
                echo_info(
                    f'More than one calculations are still running: {", ".join(str(process.pk) for process in valid)}.'
                    f' Selected {valid[index].pk}'
                )
                calcjob = valid[index]
            else:
                calcjob = valid[0]
        kwargs['calcjob'] = calcjob
        kwargs['index'] = index

        func_ = click.option(
            '--index',
            '-i',
            help='The index of the calculation to cat if multiple calculations are still running',
            default=0,
        )(func)
        return func_(*args, **kwargs)

    return wrapper


@tools.command('remotecat')
@PROCESS('calcjob')
@click.argument('fname')
@click.option('--save-to', '-s', help='Name of the file to save to')
@select_calcjob_from_work
def remotecat(calcjob: orm.CalcJobNode, fname: str, save_to: str | None) -> None:
    """
    Print the conetent of a remote file to STDOUT

    This command for printing the content of a remote file to STDOUT.
    Useful for analysing running calculations.
    """

    rfolder = calcjob.outputs.remote_folder
    if save_to is None:
        fd, temppath = tempfile.mkstemp()
    else:
        temppath = save_to

    rfolder.getfile(fname, temppath)

    with open(temppath, 'rb') as fhandle:
        copyfileobj(fhandle, sys.stdout.buffer)

    if save_to is None:
        os.close(fd)
        os.remove(temppath)


@tools.command('remotepull')
@PROCESS('calcjob')
@click.argument('dest')
@click.option(
    '--max-size',
    '-m',
    help='Maximum size of the files to be retrieved - this is passed to rsync',
)
@select_calcjob_from_work
def remotepull(calcjob: orm.CalcJobNode, dest: str, max_size: str | None) -> None:
    """
    Pull a calculation folder from the remote

    This command for pull a calculation folder to a local folder.
    `rsync` is used for doing the heavy lifting.
    """

    rfolder = calcjob.outputs.remote_folder
    cmd_args = ['rsync', '-av']

    if max_size:
        cmd_args.extend(['--max-size', max_size])

    cmd_args.append(f'{rfolder.computer.hostname}:{rfolder.get_remote_path()}/')
    if not dest.endswith('/'):
        dest = dest + '/'
    cmd_args.append(dest)

    echo_info('Running commands: {}'.format(' '.join(cmd_args)))

    completed = subprocess.run(cmd_args, check=False)
    if completed.returncode != 0:
        echo_error('Failled to pull data using rsync')
    else:
        echo_success(f'Remote folder pulled to {dest}')


@tools.command('remotetail')
@CALCULATION('calcjob')
@click.argument('fname')
def remotetail(calcjob: orm.CalcJobNode, fname: str) -> None:
    """
    Follow a file on the remote computer

    This command will launch a ssh session dedicated for following a file
    using the `tail -f` command
    """

    try:
        transport = calcjob.get_transport()
    except NotExistent as exception:
        echo_critical(repr(exception))

    remote_workdir = calcjob.get_remote_workdir()

    if not remote_workdir:
        echo_critical('no remote work directory for this calcjob, maybe the daemon did not submit it yet')

    command = tailf_command(transport, remote_workdir, fname)
    os.system(command)


@tools.command('relaxcat')
@WORKFLOW('workflow')
@click.argument('fname')
def relaxcat(workflow: orm.WorkChainNode, fname: str) -> None:
    """Cat the output of the last calculation of a finished workflow"""

    q = QueryBuilder()
    q.append(WorkChainNode, filters={'id': workflow.id})
    q.append(WorkChainNode)
    q.append(CalcJobNode, tag='calc', project=['*', 'ctime'])
    q.order_by({'calc': {'ctime': 'desc'}})
    calc, _ = q.first()

    click.Context(calcjob_outputcat).invoke(calcjob_outputcat, calcjob=calc, path=fname)


def tailf_command(transport: Any, remotedir: str, fname: str) -> str:
    """
    Specific gotocomputer string to connect to a given remote computer via
    ssh and directly go to the calculation folder and then do tail -f of the target file.
    """

    further_params = []
    if 'username' in transport._connect_args:
        further_params.append('-l {}'.format(escape_for_bash(transport._connect_args['username'])))

    if transport._connect_args.get('port'):
        further_params.append('-p {}'.format(transport._connect_args['port']))

    if transport._connect_args.get('key_filename'):
        further_params.append('-i {}'.format(escape_for_bash(transport._connect_args['key_filename'])))

    further_params_str = ' '.join(further_params)

    connect_string = (
        """ "if [ -d {escaped_remotedir} ] ;"""
        """ then cd {escaped_remotedir} ; {bash_command} -c 'tail -f {fname}' ; else echo '  ** The directory' ; """
        """echo '  ** {remotedir}' ; echo '  ** seems to have been deleted, I logout...' ; fi" """.format(
            bash_command=transport._bash_command_str,
            escaped_remotedir=f"'{remotedir}'",
            remotedir=remotedir,
            fname=fname,
        )
    )

    cmd = 'ssh -t {machine} {further_params} {connect_string}'.format(
        further_params=further_params_str,
        machine=transport._machine,
        connect_string=connect_string,
    )
    return cmd
