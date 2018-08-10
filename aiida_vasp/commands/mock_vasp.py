"""Separate cli interface for commands useful in development and testing."""
import os
import re
import click
from py import path as py_path  # pylint: disable=no-member,no-name-in-module

from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.incar import IncarParser
from aiida_vasp.io.potcar import PotcarIo
from aiida_vasp.io.poscar import PoscarParser
from aiida_vasp.io.kpoints import KpParser


def output_file(*args):
    return py_path.local(data_path(*args))


@click.command('mock-vasp')
def mock_vasp():
    """Verify input files are parseable and copy in output files."""
    from aiida.common.setup import AIIDA_CONFIG_FOLDER
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
    assert KpParser(file_path=kpoints.strpath), 'KPOINTS could not be parsed.'

    system = incar_parser.get_quantity('incar', {})['incar'].get_dict().get('system', '')
    test_case = re.findall(r'test-case:(.*)$', system)
    print("CASE", test_case)

    if not test_case:
        output_file('outcar', 'OUTCAR').copy(pwd.join('OUTCAR'))
        output_file('vasprun', 'vasprun.xml').copy(pwd.join('vasprun.xml'))
        output_file('chgcar', 'CHGCAR').copy(pwd.join('CHGCAR'))
        output_file('wavecar', 'WAVECAR').copy(pwd.join('WAVECAR'))
        output_file('eigenval', 'EIGENVAL').copy(pwd.join('EIGENVAL'))
        output_file('doscar', 'DOSCAR').copy(pwd.join('DOSCAR'))
        poscar.copy(pwd.join('CONTCAR'))
    else:
        test_case = test_case[0]
        test_data_path = py_path.local(data_path(test_case)).join('out')
        print("PATH", test_data_path)
        for out_file in test_data_path.listdir():
            out_file.copy(pwd)
