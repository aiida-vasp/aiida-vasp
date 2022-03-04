"""
Mock vasp command.

------------------
Separate cli interface for commands useful in development and testing.
"""
import os
import shutil
from pathlib import Path
import click

from aiida.cmdline.utils.decorators import with_dbenv
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.parsers.content_parsers.incar import IncarParser
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.parsers.content_parsers.kpoints import KpointsParser

from aiida_vasp.utils.mock_code import MockRegistry, MockVasp


def output_object(*args):
    return Path(data_path(*args))


@click.command('mock-vasp')
@with_dbenv()
def mock_vasp():
    """Original version of mock-vasp"""
    return _mock_vasp(False)


@click.command('mock-vasp-strict')
@with_dbenv()
def mock_vasp_strict():
    """A stricter version of mock-vasp does not allow default matching"""
    return _mock_vasp(True)


def _mock_vasp(strict_match):  # pylint: disable=too-many-statements
    """Verify input objects are parseable and copy in output objects."""
    from aiida.manage.configuration.settings import AIIDA_CONFIG_FOLDER  # pylint: disable=import-outside-toplevel
    pwd = Path().absolute()
    with open('/tmp/pung', 'w') as handler:
        handler.write('HERE')
        handler.write('pwd:' + str(pwd) + '!')
    aiida_path = Path(AIIDA_CONFIG_FOLDER)
    aiida_cfg = aiida_path / 'config.json'
    click.echo('DEBUG: AIIDA_PATH = {}'.format(os.environ.get('AIIDA_PATH')))
    click.echo('DEBUG: AIIDA_CONFIG_FOLDER = {}'.format(str(aiida_path)))
    assert aiida_path.exists()
    assert aiida_cfg.is_file()
    click.echo(aiida_cfg.read_text())
    incar = pwd / 'INCAR'
    assert incar.is_file(), 'INCAR input was not found.'

    potcar = pwd / 'POTCAR'
    assert potcar.is_file(), 'POTCAR input not found.'

    poscar = pwd / 'POSCAR'
    assert poscar.is_file(), 'POSCAR input not found.'

    kpoints = pwd / 'KPOINTS'
    assert kpoints.is_file(), 'KPOINTS input not found.'

    # Check that the input files can be parsed (as close to a validity check we can get)
    incar_parser = False
    system = ''
    with open(str(incar), 'r') as handler:
        incar_parser = IncarParser(handler=handler, validate_tags=False)
        system = incar_parser.incar.get('system', '')
    assert incar_parser, 'INCAR could not be parsed.'

    poscar_parser = False
    with open(str(poscar), 'r') as handler:
        poscar_parser = PoscarParser(handler=handler)
    assert poscar_parser, 'POSCAR could not be parsed.'

    kpoints_parser = False
    with open(str(kpoints), 'r') as handler:
        kpoints_parser = KpointsParser(handler=handler)
    assert kpoints_parser, 'KPOINTS could not be parsed.'

    #assert PotcarIo(path=str(potcar)), 'POTCAR could not be parsed.'

    try:
        test_case = system.strip().split(':')[1].strip()
    except IndexError:
        test_case = ''

    if not test_case:
        # If no test case is defined, we first try the hash-based mock registry
        mock_registry_path = os.environ.get('MOCK_CODE_BASE', data_path('.'))
        mock = MockVasp(pwd, MockRegistry(mock_registry_path))
        if mock.is_runnable:
            mock.run()
        else:
            if not strict_match:
                # Then this is a simple case - assemble the outputs from folders
                shutil.copy(output_object('outcar', 'OUTCAR'), pwd / 'OUTCAR')
                shutil.copy(output_object('vasprun', 'vasprun.xml'), pwd / 'vasprun.xml')
                shutil.copy(output_object('chgcar', 'CHGCAR'), pwd / 'CHGCAR')
                shutil.copy(output_object('wavecar', 'WAVECAR'), pwd / 'WAVECAR')
                shutil.copy(output_object('eigenval', 'EIGENVAL'), pwd / 'EIGENVAL')
                shutil.copy(output_object('doscar', 'DOSCAR'), pwd / 'DOSCAR')
                shutil.copy(output_object('basic_run', 'vasp_output'), pwd / 'vasp_output')
                shutil.copy(poscar, pwd / 'CONTCAR')
            else:
                click.echo('No matching results found but strict matching is reqested. The mock code cannot be run.')
    else:
        test_data_path = data_path(test_case, 'out')
        for out_object in Path(test_data_path).iterdir():
            shutil.copy(out_object, pwd)
