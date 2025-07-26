# pylint: disable=too-many-function-args
"""
Mock vasp command.

Separate cli interface for commands useful in development and testing.
"""

import os
import pathlib
import shutil
from typing import List

import click

from aiida_vasp.parsers.content_parsers.incar import IncarParser
from aiida_vasp.parsers.content_parsers.kpoints import KpointsParser
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.utils.mock_code import MockVasp, VaspMockRegistry, data_path


@click.command('mock-vasp-loose')
def mock_vasp_loose() -> None:
    """
    Loose version of mock-vasp that has default test data - only useful for testing and development.
    """
    return _mock_vasp(False)


@click.command('mock-vasp')
def mock_vasp() -> None:
    """
    If `MOCK_VASP_VASP_CMD` is set in the environment, it will use that command to run VASP if needed and add the
    calculation to the registry.
    """
    return _mock_vasp(True)


def _mock_vasp(strict_match: bool) -> None:  # pylint: disable=too-many-statements, too-many-locals, too-many-branches
    """
    Verify input objects are parsable and copy in output objects.
    """
    pwd = pathlib.Path().absolute()
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
        mock_registry_path = os.environ.get('MOCK_VASP_REG_BASE', data_path('.'))
        mock_registry = VaspMockRegistry(mock_registry_path)
        vasp_mock_output.append(f'MOCK PREPEND: registry search paths: {mock_registry.search_paths}\n')

        # Setup the mock code
        mock = MockVasp(pwd, mock_registry, vasp_cmd=os.environ.get('MOCK_VASP_VASP_CMD'))
        if mock.is_runnable:
            mock.run()
            detected_path = mock.registry.get_path_by_hash(mock_registry.compute_hash(pwd))
            vasp_mock_output.append(
                f'MOCK PREPEND: Using test data in path {detected_path} based detection from inputs.\n'
            )
        else:
            vasp_mock_output.append(
                'MOCK PREPEND: Using default test data in the respective folders named similar to the file name.\n'
            )
            if not strict_match:
                # Then this is a simple case - assemble the outputs from folders
                shutil.copy(data_path('outcar', 'OUTCAR'), pwd / 'OUTCAR')
                shutil.copy(data_path('vasprun', 'vasprun.xml'), pwd / 'vasprun.xml')
                shutil.copy(data_path('chgcar', 'CHGCAR'), pwd / 'CHGCAR')
                shutil.copy(data_path('wavecar', 'WAVECAR'), pwd / 'WAVECAR')
                shutil.copy(data_path('eigenval', 'EIGENVAL'), pwd / 'EIGENVAL')
                shutil.copy(data_path('doscar', 'DOSCAR'), pwd / 'DOSCAR')
                shutil.copy(data_path('basic_run', 'vasp_output'), pwd / 'vasp_output')
                shutil.copy(poscar, pwd / 'CONTCAR')
            else:
                vasp_mock_output.append(
                    'MOCK PREPEND: Caller demanded to only locate test data by input, but no match was found.\n'
                )
                stop_and_return(vasp_mock_output)
    else:
        vasp_mock_output.append('MOCK PREPEND: Using test data from folder: ' + test_case + '\n')
        test_data_path = data_path(test_case, 'out')
        for out_object in pathlib.Path(test_data_path).iterdir():
            shutil.copy(out_object, pwd)

    # Read original vasp_output as we will append mock messages to it
    vasp_output_content = []
    if vasp_output_file.exists():
        with open(vasp_output_file, 'r', encoding='utf8') as handler:
            vasp_output_content = handler.readlines()

    vasp_mock_output.append('MOCK PREPEND: Mock folder contains the following files: ' + str(os.listdir(pwd)) + '\n')
    vasp_mock_output.append('MOCK PREPEND: END ----------------------\n')
    vasp_mock_output.append('Existing VASP stdout/stderr follows:\n')

    # Make sure we add the mock details in case we need to inspect later
    with open(vasp_output_file, 'w', encoding='utf8') as handler:
        handler.write(''.join(vasp_mock_output + vasp_output_content))


def stop_and_return(vasp_mock_output: List[str]) -> None:
    """Halts mock-vasp, rebuilds the vasp_output and returns."""
    # Assemble the
    print(''.join(vasp_mock_output))
    raise RuntimeError('The mock-vasp code could not perform a clean run.')
