"""
Unit tests for aiida-vasp dryrun command.
"""

import shutil
from pathlib import Path
from unittest.mock import Mock, patch

import pytest
import yaml
from aiida.common import AttributeDict
from click.testing import CliRunner

from aiida_vasp.commands.dryrun_vasp import cmd_dryrun_vasp, dryrun_vasp, parse_ibzkpt, parse_outcar


@pytest.fixture
def cmd_params(tmp_path):
    """Common building blocks for dryrun command calls."""
    params = AttributeDict()
    params.TMP_PATH = Path(tmp_path)
    params.INPUT_DIR = Path(tmp_path / 'vasp_input')
    params.WORK_DIR = Path(tmp_path / 'vasp_work')

    # Copy existing calculation input and outputs
    shutil.copytree(Path(__file__).parent.parent / 'test_data' / 'import_calc', params.INPUT_DIR)

    return params


def run_cmd(args=None, **kwargs):
    """Run dryrun-vasp [args]."""
    runner = CliRunner()
    params = args or []
    return runner.invoke(cmd_dryrun_vasp, params, **kwargs)


def test_dryrun_help():
    """Test that help command works."""
    result = run_cmd(args=['--help'])
    assert result.exit_code == 0
    assert 'A simple tool to dryrun a VASP calculation' in result.output


def test_dryrun_default_parameters():
    """Test dryrun with default parameters."""
    with patch('aiida_vasp.commands.dryrun_vasp.dryrun_vasp') as mock_dryrun:
        mock_dryrun.return_value = {'num_kpoints': 10, 'num_bands': 20}
        result = run_cmd(args=[])
        assert result.exit_code == 0
        mock_dryrun.assert_called_once()


def test_dryrun_custom_input_dir(cmd_params):
    """Test dryrun with custom input directory."""
    with patch('aiida_vasp.commands.dryrun_vasp.dryrun_vasp') as mock_dryrun:
        mock_dryrun.return_value = {'num_kpoints': 10, 'num_bands': 20}

        result = run_cmd(args=['--input-dir', cmd_params.INPUT_DIR])
        assert result.exit_code == 0
        mock_dryrun.assert_called_once_with(
            input_dir=cmd_params.INPUT_DIR, vasp_exe='vasp_std', timeout=10, work_dir=None, keep=False, force=False
        )


def test_dryrun_custom_vasp_exe(cmd_params):
    """Test dryrun with custom VASP executable."""
    with patch('aiida_vasp.commands.dryrun_vasp.dryrun_vasp') as mock_dryrun:
        mock_dryrun.return_value = {'num_kpoints': 10, 'num_bands': 20}

        result = run_cmd(args=['--input-dir', cmd_params.INPUT_DIR, '--vasp-exe', 'custom_vasp'])
        assert result.exit_code == 0
        mock_dryrun.assert_called_once_with(
            input_dir=cmd_params.INPUT_DIR, vasp_exe='custom_vasp', timeout=10, work_dir=None, keep=False, force=False
        )


def test_dryrun_custom_timeout(cmd_params):
    """Test dryrun with custom timeout."""
    with patch('aiida_vasp.commands.dryrun_vasp.dryrun_vasp') as mock_dryrun:
        mock_dryrun.return_value = {'num_kpoints': 10, 'num_bands': 20}

        result = run_cmd(args=['--input-dir', cmd_params.INPUT_DIR, '--timeout', '30'])
        assert result.exit_code == 0
        mock_dryrun.assert_called_once_with(
            input_dir=cmd_params.INPUT_DIR, vasp_exe='vasp_std', timeout=30, work_dir=None, keep=False, force=False
        )


def test_dryrun_custom_work_dir(cmd_params):
    """Test dryrun with custom work directory."""
    with patch('aiida_vasp.commands.dryrun_vasp.dryrun_vasp') as mock_dryrun:
        mock_dryrun.return_value = {'num_kpoints': 10, 'num_bands': 20}

        result = run_cmd(args=['--input-dir', cmd_params.INPUT_DIR, '--work-dir', cmd_params.WORK_DIR])
        assert result.exit_code == 0
        mock_dryrun.assert_called_once_with(
            input_dir=cmd_params.INPUT_DIR,
            vasp_exe='vasp_std',
            timeout=10,
            work_dir=str(cmd_params.WORK_DIR),
            keep=False,
            force=False,
        )


def test_dryrun_keep_files(cmd_params):
    """Test dryrun with keep flag."""
    with patch('aiida_vasp.commands.dryrun_vasp.dryrun_vasp') as mock_dryrun:
        mock_dryrun.return_value = {'num_kpoints': 10, 'num_bands': 20}

        result = run_cmd(args=['--input-dir', cmd_params.INPUT_DIR, '--keep'])
        assert result.exit_code == 0
        mock_dryrun.assert_called_once_with(
            input_dir=cmd_params.INPUT_DIR, vasp_exe='vasp_std', timeout=10, work_dir=None, keep=True, force=False
        )


def test_dryrun_force_flag(cmd_params):
    """Test dryrun with force flag."""
    with patch('aiida_vasp.commands.dryrun_vasp.dryrun_vasp') as mock_dryrun:
        mock_dryrun.return_value = {'num_kpoints': 10, 'num_bands': 20}

        result = run_cmd(args=['--input-dir', cmd_params.INPUT_DIR, '--force'])
        assert result.exit_code == 0
        mock_dryrun.assert_called_once_with(
            input_dir=cmd_params.INPUT_DIR, vasp_exe='vasp_std', timeout=10, work_dir=None, keep=False, force=True
        )


def test_dryrun_creates_yaml_output(cmd_params):
    """Test that dryrun creates yaml output file."""
    with patch('aiida_vasp.commands.dryrun_vasp.dryrun_vasp') as mock_dryrun:
        mock_result = {'num_kpoints': 10, 'num_bands': 20, 'NGX': 24}
        mock_dryrun.return_value = mock_result

        result = run_cmd(args=['--input-dir', cmd_params.INPUT_DIR])
        assert result.exit_code == 0

        # Check that dryrun.yaml was created
        yaml_file = Path(cmd_params.INPUT_DIR) / 'dryrun.yaml'
        assert yaml_file.exists()

        # Check content
        with open(yaml_file) as f:
            data = yaml.safe_load(f)
        assert data == mock_result


def test_dryrun_vasp_same_directories():
    """Test that dryrun_vasp fails when input and work directories are the same."""
    with pytest.raises(ValueError, match='The working directory cannot be the input directory'):
        dryrun_vasp(input_dir='/test/path', work_dir='/test/path')


def test_dryrun_vasp_existing_work_dir_no_force(tmp_path):
    """Test that dryrun_vasp fails when work directory exists and force=False."""
    input_dir = tmp_path / 'input'
    work_dir = tmp_path / 'work'
    input_dir.mkdir()
    work_dir.mkdir()  # Create existing work directory

    with pytest.raises(FileExistsError, match=r'Working directory .* exists already'):
        dryrun_vasp(input_dir=str(input_dir), work_dir=str(work_dir), force=False)


def test_dryrun_vasp_existing_work_dir_with_force(cmd_params):
    """Test that dryrun_vasp succeeds when work directory exists and force=True."""
    with (
        patch('subprocess.Popen') as mock_popen,
        patch('time.sleep'),
        patch('aiida_vasp.commands.dryrun_vasp.parse_outcar') as mock_parse_outcar,
        patch('aiida_vasp.commands.dryrun_vasp.parse_ibzkpt') as mock_parse_ibzkpt,
    ):
        # Mock process
        mock_process = Mock()
        mock_process.poll.return_value = None
        mock_popen.return_value = mock_process

        mock_parse_outcar.return_value = {'num_kpoints': 10}
        mock_parse_ibzkpt.return_value = []

        # Should not raise an error
        result = dryrun_vasp(
            input_dir=str(cmd_params.INPUT_DIR),
            work_dir=str(cmd_params.WORK_DIR),
            force=True,
            timeout=1,  # Short timeout for testing
        )
        assert isinstance(result, dict)


def test_parse_outcar(cmd_params):
    """Test OUTCAR parsing functionality."""

    outcar_file = cmd_params.INPUT_DIR / 'OUTCAR'
    result = parse_outcar(outcar_file)

    assert result['POTCARS'] == ['PAW_PBE Si 05Jan2001']
    assert result['num_kpoints'] == 20
    assert result['num_bands'] == 9
    assert result['NGX'] == 16
    assert result['NGY'] == 16
    assert result['NGZ'] == 16
    assert result['num_plane_waves'] == 4096
    assert result['plane_waves_min_max'] == [359.0, 338.0]
    assert result['max_ram_rank0'] == 36636.0
    assert result['mem_base'] == 30000.0
    assert 'kpoints_and_weights' in result


def test_parse_ibzkpt(cmd_params):
    """Test IBZKPT parsing functionality."""
    ibzkpt_file = cmd_params.INPUT_DIR / 'IBZKPT'
    result = parse_ibzkpt(ibzkpt_file)

    assert len(result) == 20
    # Check that weights are normalized
    total_weight = sum(point[3] for point in result)
    assert abs(total_weight - 1.0) < 1e-10

    # Check first point
    assert result[0][:3] == [0.0, 0.0, 0.0]


def test_parse_outcar_empty_file(tmp_path):
    """Test OUTCAR parsing with empty file."""
    outcar_file = tmp_path / 'OUTCAR'
    outcar_file.write_text('')

    result = parse_outcar(outcar_file)

    assert result['POTCARS'] == []
    assert 'num_kpoints' not in result
