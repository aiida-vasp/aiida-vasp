"""
Mock vasp command.

------------------
Separate cli interface for commands useful in development and testing.
"""
import os
import click
from py import path as py_path  # pylint: disable=no-member,no-name-in-module

from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.parsers.file_parsers.incar import IncarParser
from aiida_vasp.parsers.file_parsers.potcar import PotcarIo
from aiida_vasp.parsers.file_parsers.poscar import PoscarParser
from aiida_vasp.parsers.file_parsers.kpoints import KpointsParser


def output_file(*args):
    return py_path.local(data_path(*args))


@click.command('mock-vasp')
def mock_vasp():
    """Verify input files are parseable and copy in output files."""
    from aiida.manage.configuration.settings import AIIDA_CONFIG_FOLDER  # pylint: disable=import-outside-toplevel

    pwd = py_path.local('.')
    aiida_path = py_path.local(AIIDA_CONFIG_FOLDER)
    aiida_cfg = aiida_path.join('config.json')
    click.echo('DEBUG: AIIDA_PATH = {}'.format(os.environ.get('AIIDA_PATH')))
    click.echo('DEBUG: AIIDA_CONFIG_FOLDER = {}'.format(aiida_path.strpath))
    assert aiida_path.isdir()
    assert aiida_cfg.isfile()
    click.echo(aiida_cfg.read())
    incar = pwd.join('INCAR')
    assert incar.isfile(), 'INCAR input file was not found.'

    potcar = pwd.join('POTCAR')
    assert potcar.isfile(), 'POTCAR input file not found.'

    poscar = pwd.join('POSCAR')
    assert poscar.isfile(), 'POSCAR input file not found.'

    kpoints = pwd.join('KPOINTS')
    assert kpoints.isfile(), 'KPOINTS input file not found.'

    incar_parser = IncarParser(file_path=incar.strpath)
    assert incar_parser, 'INCAR could not be parsed.'
    assert PotcarIo(path=potcar.strpath), 'POTCAR could not be parsed.'
    assert PoscarParser(file_path=poscar.strpath), 'POSCAR could not be parsed.'
    assert KpointsParser(file_path=kpoints.strpath), 'KPOINTS could not be parsed.'

    system = incar_parser.incar.get('system', '')
    try:
        test_case = system.strip().split(':')[1].strip()
    except IndexError:
        test_case = ''
    if not test_case:
        output_file('outcar', 'OUTCAR').copy(pwd.join('OUTCAR'))
        output_file('vasprun', 'vasprun.xml').copy(pwd.join('vasprun.xml'))
        output_file('chgcar', 'CHGCAR').copy(pwd.join('CHGCAR'))
        output_file('wavecar', 'WAVECAR').copy(pwd.join('WAVECAR'))
        output_file('eigenval', 'EIGENVAL').copy(pwd.join('EIGENVAL'))
        output_file('doscar', 'DOSCAR').copy(pwd.join('DOSCAR'))
        poscar.copy(pwd.join('CONTCAR'))
    else:
        test_data_path = data_path(test_case, 'out')
        for out_file in py_path.local(test_data_path).listdir():
            out_file.copy(pwd)
