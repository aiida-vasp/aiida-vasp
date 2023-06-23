"""
Mock vasp command.

------------------
Separate cli interface for commands useful in development and testing.
"""
import os
import shutil
from pathlib import Path
import click

from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.parsers.file_parsers.incar import IncarParser
from aiida_vasp.parsers.file_parsers.potcar import PotcarIo
from aiida_vasp.parsers.file_parsers.poscar import PoscarParser
from aiida_vasp.parsers.file_parsers.kpoints import KpointsParser


def output_file(*args):
    return Path(data_path(*args))


@click.command('mock-vasp')
def mock_vasp():
    """Verify input files are parseable and copy in output files."""
    from aiida.manage.configuration.settings import AIIDA_CONFIG_FOLDER  # pylint: disable=import-outside-toplevel
    pwd = Path().absolute()
    aiida_path = Path(AIIDA_CONFIG_FOLDER)
    aiida_cfg = aiida_path / 'config.json'
    click.echo('DEBUG: AIIDA_PATH = {}'.format(os.environ.get('AIIDA_PATH')))
    click.echo('DEBUG: AIIDA_CONFIG_FOLDER = {}'.format(str(aiida_path)))
    assert aiida_path.exists()
    assert aiida_cfg.is_file()
    click.echo(aiida_cfg.read_text())
    incar = pwd / 'INCAR'
    assert incar.is_file(), 'INCAR input file was not found.'

    potcar = pwd / 'POTCAR'
    assert potcar.is_file(), 'POTCAR input file not found.'

    poscar = pwd / 'POSCAR'
    assert poscar.is_file(), 'POSCAR input file not found.'

    kpoints = pwd / 'KPOINTS'
    assert kpoints.is_file(), 'KPOINTS input file not found.'
    incar_parser = IncarParser(file_path=str(incar))
    assert incar_parser, 'INCAR could not be parsed.'
    assert PotcarIo(path=str(potcar)), 'POTCAR could not be parsed.'
    assert PoscarParser(file_path=str(poscar)), 'POSCAR could not be parsed.'
    assert KpointsParser(file_path=str(kpoints)), 'KPOINTS could not be parsed.'

    system = incar_parser.incar.get('system', '')
    try:
        test_case = system.strip().split(':')[1].strip()
    except IndexError:
        test_case = ''
    if not test_case:
        shutil.copy(output_file('outcar', 'OUTCAR'), pwd / 'OUTCAR')
        shutil.copy(output_file('vasprun', 'vasprun.xml'), pwd / 'vasprun.xml')
        shutil.copy(output_file('chgcar', 'CHGCAR'), pwd / 'CHGCAR')
        shutil.copy(output_file('wavecar', 'WAVECAR'), pwd / 'WAVECAR')
        shutil.copy(output_file('eigenval', 'EIGENVAL'), pwd / 'EIGENVAL')
        shutil.copy(output_file('doscar', 'DOSCAR'), pwd / 'DOSCAR')
        shutil.copy(output_file('basic_run', 'vasp_output'), pwd / 'vasp_output')
        shutil.copy(poscar, pwd / 'CONTCAR')
    else:
        test_data_path = data_path(test_case, 'out')
        for out_file in Path(test_data_path).iterdir():
            shutil.copy(out_file, pwd)
