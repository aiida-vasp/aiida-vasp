"""Separate cli interface for commands useful in development and testing."""
import click
from py import path as py_path  # pylint: disable=no-member,no-name-in-module

from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.incar import IncarIo
from aiida_vasp.io.potcar import PotcarIo
from aiida_vasp.io.poscar import PoscarParser
from aiida_vasp.io.kpoints import KpParser


def output_file(*args):
    return py_path.local(data_path(*args))


@click.command('mock-vasp')
def mock_vasp():
    """Verify input files are parseable and copy in output files."""
    pwd = py_path.local('.')

    incar = pwd.join('INCAR')
    assert incar.isfile(), 'INCAR input file was not found.'

    potcar = pwd.join('POTCAR')
    assert potcar.isfile(), 'POTCAR input file not found.'

    poscar = pwd.join('POSCAR')
    assert poscar.isfile(), 'POSCAR input file not found.'

    kpoints = pwd.join('KPOINTS')
    assert kpoints.isfile(), 'KPOINTS input file not found.'

    assert IncarIo(file_path=incar.strpath), 'INCAR could not be parsed.'
    assert PotcarIo(path=potcar.strpath), 'POTCAR could not be parsed.'
    assert PoscarParser(file_path=poscar.strpath), 'POSCAR could not be parsed.'
    assert KpParser(file_path=kpoints.strpath), 'KPOINTS could not be parsed.'

    output_file('outcar', 'OUTCAR').copy(pwd.join('OUTCAR'))
    output_file('vasprun', 'vasprun.xml').copy(pwd.join('vasprun.xml'))
    output_file('chgcar', 'CHGCAR').copy(pwd.join('CHGCAR'))
    output_file('wavecar', 'WAVECAR').copy(pwd.join('WAVECAR'))
    output_file('eigenval', 'EIGENVAL').copy(pwd.join('EIGENVAL'))
    output_file('doscar', 'DOSCAR').copy(pwd.join('DOSCAR'))
    poscar.copy(pwd.join('CONTCAR'))
