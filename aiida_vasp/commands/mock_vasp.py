# pylint: disable=too-many-function-args
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
from aiida_vasp.parsers.content_parsers.incar import IncarParser
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.parsers.content_parsers.kpoints import KpointsParser

from aiida_vasp.utils.mock_code import VaspMockRegistry, MockVasp


def output_object(*args):
    return Path(data_path(*args))


@click.command('mock-vasp')
def mock_vasp():
    """Original version of mock-vasp"""
    return _mock_vasp(False)


@click.command('mock-vasp-strict')
def mock_vasp_strict():
    """A stricter version of mock-vasp does not allow default matching"""
    return _mock_vasp(True)


def _mock_vasp(strict_match):  # pylint: disable=too-many-statements, too-many-locals, too-many-branches
    """Verify input objects are parseable and copy in output objects."""
    pwd = Path().absolute()
    vasp_mock_output = []
    vasp_output_file = pwd / 'vasp_output'
    vasp_mock_output.append('MOCK PREPEND: START ----------------------\n')
    vasp_mock_output.append('MOCK PREPEND: Mock directory: ' + str(pwd) + '\n')

    incar = pwd / 'INCAR'
    if not incar.is_file():
        vasp_mock_output.append('MOCK PREPEND: INCAR not found.\n')
        stop_and_return(vasp_mock_output)

    potcar = pwd / 'POTCAR'
    if not potcar.is_file():
        vasp_mock_output.append('MOCK PREPEND: POTCAR not found.\n')
        stop_and_return(vasp_mock_output)

    poscar = pwd / 'POSCAR'
    if not poscar.is_file():
        vasp_mock_output.append('MOCK PREPEND: POSCAR not found.\n')
        stop_and_return(vasp_mock_output)

    kpoints = pwd / 'KPOINTS'
    if not kpoints.is_file():
        vasp_mock_output.append('MOCK PREPEND: KPOINTS not found.\n')
        stop_and_return(vasp_mock_output)

    # Check that the input files can be parsed (as close to a validity check we can get)
    incar_parser = False
    system = ''
    with open(str(incar), 'r', encoding='utf8') as handler:
        incar_parser = IncarParser(handler=handler, validate_tags=False)
        system = incar_parser.incar.get('system', '')
    if not incar_parser:
        vasp_mock_output.append('MOCK PREPEND: INCAR could not be parsed.\n')
        stop_and_return(vasp_mock_output)

    poscar_parser = False
    with open(str(poscar), 'r', encoding='utf8') as handler:
        poscar_parser = PoscarParser(handler=handler)
    if not poscar_parser:
        vasp_mock_output.append('MOCK PREPEND: POSCAR could not be parsed.\n')
        stop_and_return(vasp_mock_output)

    kpoints_parser = False
    with open(str(kpoints), 'r', encoding='utf8') as handler:
        kpoints_parser = KpointsParser(handler=handler)
    if not kpoints_parser:
        vasp_mock_output.append('MOCK PREPEND: KPOINTS could not be parsed.\n')
        stop_and_return(vasp_mock_output)

    try:
        test_case = system.strip().split(':')[1].strip()
    except IndexError:
        test_case = ''

    if not test_case:
        vasp_mock_output.append('MOCK PREPEND: Trying to detect test case using registry or reverting to default.\n')
        # If no test case is defined, we first try the hash-based mock registry
        mock_registry_path = os.environ.get('VASP_MOCK_CODE_BASE', data_path('.'))
        mock_registry = VaspMockRegistry(mock_registry_path)
        vasp_mock_output.append(f'MOCK PREPEND: registry search paths: {mock_registry.search_paths}\n')
        mock = MockVasp(pwd, mock_registry)
        if mock.is_runnable:
            detected_path = mock.registry.get_path_by_hash(mock_registry.compute_hash(pwd))
            vasp_mock_output.append(f'MOCK PREPEND: Using test data in path {detected_path} based detection from inputs.\n')
            mock.run()
        else:
            vasp_mock_output.append('MOCK PREPEND: Using default test data in the respective folders named similar to the file name.\n')
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
                vasp_mock_output.append('MOCK PREPEND: Caller demanded to only locate test data by input, but no match was found.\n')
                stop_and_return(vasp_mock_output)
    else:
        vasp_mock_output.append('MOCK PREPEND: Using test data from folder: ' + test_case + '\n')
        test_data_path = data_path(test_case, 'out')
        for out_object in Path(test_data_path).iterdir():
            shutil.copy(out_object, pwd)

    # Read original vasp_output as we will append mock messages to it
    if vasp_output_file.exists():
        with open(vasp_output_file, 'r', encoding='utf8') as handler:
            vasp_output = handler.readlines()

    vasp_mock_output.append('MOCK PREPEND: Mock folder contains the following files: ' + str(os.listdir(pwd)) + '\n')
    vasp_mock_output.append('MOCK PREPEND: END ----------------------\n')
    vasp_mock_output.append('Existing VASP stdout/stderr follows:\n')

    # Make sure we add the mock details in case we need to inspect later
    with open(vasp_output_file, 'w', encoding='utf8') as handler:
        handler.write(''.join(vasp_mock_output + vasp_output))


def stop_and_return(vasp_mock_output):
    """Halts mock-vasp, rebuilds the vasp_output and returns."""
    # Assemble the
    print(''.join(vasp_mock_output))
    raise RuntimeError('The mock-vasp code could not perform a clean run.')
