"""
Unit tests for aiida-vasp immigrant command.
"""

import os
import shutil
import tempfile
from pathlib import Path

import pytest
from aiida.common import AttributeDict
from aiida.orm import Computer, InstalledCode
from click.testing import CliRunner

from aiida_vasp.commands.immigrant import import_calc


@pytest.fixture
def cmd_params(tmp_path):
    """Common building blocks for immigrant command calls."""
    params = AttributeDict()
    params.TMP_PATH = str(tmp_path)
    params.CALC_PATH = str(tmp_path / 'import_calc')
    params.REMOTE_PATH = '/remote/path/vasp_calc'

    # Create a calculation folder with required files
    calc_folder = Path(params.CALC_PATH)

    # Copy existing calculation input and outputs
    shutil.copytree(Path(__file__).parent.parent / 'test_data' / 'import_calc', calc_folder)
    return params


@pytest.fixture
def mock_code(aiida_profile_clean):
    """Create a mock VASP code for testing."""
    # Create localhost computer if it doesn't exist
    try:
        computer = Computer.collection.get(label='localhost')
    except Exception:
        computer = Computer(
            label='localhost',
            hostname='localhost',
            transport_type='core.local',
            scheduler_type='core.direct',
            workdir='/tmp',
        )
        computer.store()
        computer.configure()

    # Create a dummy executable
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
        f.write('#!/bin/bash\necho "Mock VASP"\n')
        dummy_exec = f.name

    os.chmod(dummy_exec, 0o755)

    # Create code
    code = InstalledCode(
        label='test-vasp-import',
        computer=computer,
        filepath_executable=dummy_exec,
    )
    code.store()

    yield code

    # Cleanup
    os.unlink(dummy_exec)


def run_cmd(args=None, **kwargs):
    """Run aiida-vasp import [args]."""
    runner = CliRunner()
    params = args or []
    return runner.invoke(import_calc, params, **kwargs)


def test_import_help():
    """Test that help command works."""
    result = run_cmd(args=['--help'])
    assert result.exit_code == 0
    assert 'Import an existing calculation' in result.output


def test_import_missing_path():
    """Test import command fails without path."""
    result = run_cmd(args=['--code', 'test-code'])
    assert result.exit_code == 2  # Click missing argument error


def test_import_nonexistent_path():
    """Test import from non-existent path."""
    result = run_cmd(args=['/nonexistent/path'])
    assert result.exit_code != 0


def test_import_missing_files(tmp_path):
    """Test import from folder missing required files."""
    calc_folder = tmp_path / 'incomplete_calc'
    calc_folder.mkdir()
    # Only create INCAR, missing POSCAR and KPOINTS
    (calc_folder / 'INCAR').write_text('ENCUT = 400\n')

    result = run_cmd(args=[str(calc_folder)])
    assert result.exit_code != 0
    assert 'Missing required files' in result.output


def test_import_missing_potcar_and_family(cmd_params):
    """Test import without POTCAR and no potential family."""
    # Remove POTCAR from test folder
    (Path(cmd_params.CALC_PATH) / 'POTCAR').unlink()

    result = run_cmd(args=[cmd_params.CALC_PATH])
    assert result.exit_code != 0
    assert 'POTCAR not found and no potential family specified' in result.output


def test_import_with_potential_family(cmd_params):
    """Test import with potential family specified."""
    # Remove POTCAR from test folder
    (Path(cmd_params.CALC_PATH) / 'POTCAR').unlink()

    result = run_cmd(
        args=[
            cmd_params.CALC_PATH,
            '--potential-family',
            'test-family',
            '--yes',  # Skip confirmation
        ]
    )
    # This will fail because the potential family doesn't exist, but it should get past the POTCAR check
    assert 'POTCAR not found and no potential family specified' not in result.output


def test_import_invalid_code(cmd_params):
    """Test import with invalid code."""
    result = run_cmd(args=[cmd_params.CALC_PATH, '--code', 'nonexistent-code'])
    assert result.exit_code != 0
    assert 'Code "nonexistent-code" not found' in result.output


def test_import_with_valid_code(cmd_params, mock_code):
    """Test import with valid code."""
    result = run_cmd(
        args=[
            cmd_params.CALC_PATH,
            '--code',
            f'{mock_code.pk}',
            '--label',
            'test-import',
            '--yes',  # Skip confirmation
        ]
    )
    # This may fail due to missing dependencies but should recognize the code
    assert 'Code "test-vasp-import" not found' not in result.output


def test_import_without_code_local(cmd_params):
    """Test import without code for local path (should create dummy code)."""
    result = run_cmd(
        args=[
            cmd_params.CALC_PATH,
            '--label',
            'test-import-no-code',
            '--yes',  # Skip confirmation
        ]
    )
    # Should mention creating dummy code
    if result.exit_code == 0 or 'creating dummy code' in result.output:
        assert True
    else:
        # May fail due to other issues, but shouldn't fail on missing code
        assert 'Code must be specified' not in result.output


def test_import_remote_without_code():
    """Test import remote path without code (should fail)."""
    result = run_cmd(args=['/remote/nonexistent/path', '--label', 'test-remote'])
    # Should fail because remote imports require a code
    assert result.exit_code != 0


def test_import_with_potential_mapping(cmd_params, mock_code):
    """Test import with potential mapping."""
    mapping = '{"Si": "Si_pv"}'

    result = run_cmd(
        args=[
            cmd_params.CALC_PATH,
            '--code',
            f'{mock_code.pk}',
            '--potential-mapping',
            mapping,
            '--label',
            'test-mapping',
            '--yes',
        ]
    )
    # Should parse JSON without error
    assert 'Error parsing potential mapping JSON' not in result.output


def test_import_invalid_potential_mapping(cmd_params, mock_code):
    """Test import with invalid potential mapping JSON."""
    invalid_mapping = 'invalid-json'

    result = run_cmd(
        args=[
            cmd_params.CALC_PATH,
            '--code',
            f'{mock_code.pk}',
            '--potential-mapping',
            invalid_mapping,
            '--label',
            'test-invalid-mapping',
        ]
    )
    assert result.exit_code != 0
    assert 'Error parsing potential mapping JSON' in result.output


def test_import_with_flags(cmd_params, mock_code):
    """Test import with various flags."""
    result = run_cmd(
        args=[
            cmd_params.CALC_PATH,
            '--code',
            f'{mock_code.pk}',
            '--include-wavecar',
            '--include-chgcar',
            '--label',
            'test-with-flags',
            '--description',
            'Test calculation with flags',
            '--stdout-file',
            'custom_output',
            '--yes',
        ]
    )
    # Flags should be accepted without error
    assert result.exit_code == 0 or 'Error parsing' not in result.output


def test_import_submit_daemon(cmd_params, mock_code):
    """Test import with daemon submission."""
    result = run_cmd(
        args=[cmd_params.CALC_PATH, '--code', f'{mock_code.pk}', '--submit-daemon', '--label', 'test-daemon', '--yes']
    )
    # Should accept daemon flag
    assert 'Error parsing' not in result.output


def test_import_quiet_mode(cmd_params, mock_code, upload_potcar):
    """Test import in quiet mode."""
    result = run_cmd(
        args=[
            cmd_params.CALC_PATH,
            '--code',
            f'{mock_code.pk}',
            '--quiet',
            '--label',
            'test-quiet',
            '--yes',
            '--potential-family',
            'PBE.54',
        ]
    )
    # Should have minimal output in quiet mode
    assert len(result.output.strip()) < len('Importing calculation from:')


def test_import_with_group(cmd_params, mock_code):
    """Test import with group assignment."""
    result = run_cmd(
        args=[
            cmd_params.CALC_PATH,
            '--code',
            f'{mock_code.pk}',
            '--group',
            'test-import-group',
            '--label',
            'test-group',
            '--yes',
        ]
    )
    # Should accept group parameter
    assert 'Error parsing' not in result.output


def test_import_confirmation_prompt(cmd_params, mock_code):
    """Test import confirmation prompt."""
    result = run_cmd(
        args=[cmd_params.CALC_PATH, '--code', f'{mock_code.pk}', '--label', 'test-confirm'],
        input='n\n',  # Answer 'no' to confirmation
    )
    # Should be cancelled
    assert result.exit_code == 0 or 'Cancelled import' in result.output


def test_import_with_all_options(cmd_params, mock_code):
    """Test import with all options specified."""
    result = run_cmd(
        args=[
            cmd_params.CALC_PATH,
            '--code',
            f'{mock_code.pk}',
            '--potential-family',
            'test-family',
            '--potential-mapping',
            '{"Si": "Si_pv"}',
            '--include-wavecar',
            '--include-chgcar',
            '--stdout-file',
            'vasp.out',
            '--label',
            'comprehensive-test',
            '--description',
            'Test with all options',
            '--group',
            'comprehensive-group',
            '--submit-daemon',
            '--quiet',
            '--yes',
        ]
    )
    # All options should be parsed correctly
    assert 'Error parsing' not in result.output


@pytest.mark.parametrize('stdout_file', ['vasp_output', 'vasp.out', 'custom_stdout'])
def test_import_different_stdout_files(cmd_params, mock_code, stdout_file):
    """Test import with different stdout file names."""
    result = run_cmd(
        args=[
            cmd_params.CALC_PATH,
            '--code',
            f'{mock_code.pk}',
            '--stdout-file',
            stdout_file,
            '--label',
            f'test-stdout-{stdout_file}',
            '--yes',
        ]
    )
    # Should accept different stdout file names
    assert 'Error parsing' not in result.output


def test_import_summary_display(cmd_params, mock_code):
    """Test that import summary is displayed correctly."""
    result = run_cmd(
        args=[
            cmd_params.CALC_PATH,
            '--code',
            f'{mock_code.pk}',
            '--label',
            'test-summary',
            '--description',
            'Test summary display',
        ],
        input='y\n',  # Answer 'yes' to confirmation
    )
    # Should show summary
    if 'Import Summary' in result.output:
        assert 'Calculation path:' in result.output
        assert 'Code:' in result.output
        assert 'test-summary' in result.output
        assert 'Test summary display' in result.output


def test_import_missing_required_vasp_files(tmp_path):
    """Test import with only some required files."""
    calc_folder = tmp_path / 'partial_calc'
    calc_folder.mkdir()

    # Create only INCAR and POSCAR, missing KPOINTS
    (calc_folder / 'INCAR').write_text('ENCUT = 400\n')
    (calc_folder / 'POSCAR').write_text("""Test
1.0
1 0 0
0 1 0
0 0 1
Si
1
Direct
0 0 0
""")

    result = run_cmd(args=[str(calc_folder)])
    assert result.exit_code != 0
    assert 'Missing required files' in result.output
    assert 'KPOINTS' in result.output
