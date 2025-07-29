"""Unit tests for VASP calculation monitoring functions."""

import os
import time
from pathlib import Path

import pytest
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.calcs.monitors import monitor_loop_time, monitor_stdout
from aiida_vasp.calcs.vasp import VaspCalculation


class MockFileStatResult:
    """Mock object for file stat results."""

    def __init__(self, size=1024, mtime=None):
        """Initialize mock file stat result.

        :param size: File size in bytes
        :param mtime: File modification time (defaults to current time)
        """
        self.st_size = size
        self.st_mtime = mtime if mtime is not None else time.time()


class MockTransport:
    """Mock transport object for testing."""

    def __init__(self):
        """Initialize mock transport."""
        self.get_attribute_return = None
        self.exec_command_wait_return = None
        self.get_attribute_calls = []
        self.exec_command_wait_calls = []

    def reset(self):
        """Reset the status"""
        self.get_attribute_return = None
        self.exec_command_wait_return = None
        self.get_attribute_calls = []
        self.exec_command_wait_calls = []

    def get_attribute(self, filename):
        """Mock get_attribute method."""
        self.get_attribute_calls.append(filename)
        return self.get_attribute_return

    def exec_command_wait(self, command):
        """Mock exec_command_wait method."""
        self.exec_command_wait_calls.append(command)
        if self.exec_command_wait_return is None:
            return (0, '', '')
        return self.exec_command_wait_return


class MockProcessClass:
    """Mock process class with VASP output name."""

    def __init__(self, vasp_output='vasp_output'):
        """Initialize with VASP output filename."""
        self._VASP_OUTPUT = vasp_output


class MockNode:
    """Mock CalcJobNode for testing."""

    def __init__(self, walltime_limit=None, vasp_output='vasp_output'):
        """Initialize mock node."""
        self.process_class = MockProcessClass(vasp_output)
        self.inputs = AttributeDict()
        self.options = {'max_wallclock_seconds': walltime_limit}

    def get_option(self, key):
        """Mocked get_option method"""
        return self.options[key]

    def get_remote_workdir(self):
        """Mocked get_remote_workdir method"""
        return os.getcwd()  # Use current working directory for simplicity


@pytest.fixture
def mock_transport():
    """Provide a mock transport object."""
    return MockTransport()


@pytest.fixture
def mock_node():
    """Provide a mock CalcJobNode."""
    return MockNode()


class TestMonitorStdout:
    """Test cases for monitor_stdout function."""

    def test_monitor_stdout_normal_size(self, mock_node, mock_transport):
        """Test monitor_stdout with normal file size."""
        # Setup mocks
        mock_transport.get_attribute_return = MockFileStatResult(size=1024 * 1024)  # 1MB

        # Call function
        result = monitor_stdout(mock_node, mock_transport, size_threshold_mb=5)

        # Assertions
        stdout_path = str(Path(os.getcwd()) / VaspCalculation._VASP_OUTPUT)
        assert result is None
        assert mock_transport.get_attribute_calls == [stdout_path]
        assert mock_transport.exec_command_wait_calls == []

    def test_monitor_stdout_oversized_file(self, mock_node, mock_transport):
        """Test monitor_stdout with oversized file triggering truncation."""
        # Setup mocks
        file_size = 10 * 1024 * 1024  # 10MB
        mock_transport.get_attribute_return = MockFileStatResult(size=file_size)

        # Call function
        result = monitor_stdout(mock_node, mock_transport, size_threshold_mb=5)

        # Assertions
        assert result is not None
        assert 'Very large stdout detected' in result
        assert '10.00 MB' in result
        assert 'potential critical crash' in result

        stdout_path = str(Path(os.getcwd()) / VaspCalculation._VASP_OUTPUT)
        assert mock_transport.get_attribute_calls == [stdout_path]
        assert mock_transport.exec_command_wait_calls == [f'truncate -s 5M {stdout_path}']

        # Apply a custom threshold
        mock_transport.reset()
        mock_transport.get_attribute_return = MockFileStatResult(size=file_size)
        result = monitor_stdout(mock_node, mock_transport, size_threshold_mb=15)

        assert result is None
        assert not mock_transport.exec_command_wait_calls

        # Also pass if exactly at threshold
        # Apply a custom threshold
        mock_transport.reset()
        mock_transport.get_attribute_return = MockFileStatResult(size=file_size)
        result = monitor_stdout(mock_node, mock_transport, size_threshold_mb=10)
        assert result is None
        assert not mock_transport.exec_command_wait_calls


class TestMonitorLoopTime:
    """Test cases for monitor_loop_time function."""

    def test_monitor_loop_time_no_walltime_limit(self, mock_transport):
        """Test monitor_loop_time when no walltime limit is set."""
        # Setup node without walltime limit
        mock_node = MockNode(walltime_limit=None)

        # Call function
        result = monitor_loop_time(mock_node, mock_transport)

        # Assertions
        assert result is None
        assert mock_transport.exec_command_wait_calls == []

    def test_monitor_loop_time_no_loop_entries(self, mock_transport):
        """Test monitor_loop_time when no LOOP entries are found in OUTCAR."""
        # Setup node with walltime
        mock_node = MockNode(walltime_limit=3600)  # 1 hour
        mock_transport.exec_command_wait_return = (1, '', 'grep: no match')

        # Call function
        result = monitor_loop_time(mock_node, mock_transport)

        # Assertions
        assert result is None
        outcar_path = Path(os.getcwd()) / 'OURCAR'
        assert mock_transport.exec_command_wait_calls == [f"grep 'LOOP:' {outcar_path}"]

    def test_monitor_loop_time_fast_loops(self, mock_transport):
        """Test monitor_loop_time with fast electronic loops."""
        # Setup node with walltime
        mock_node = MockNode(walltime_limit=3600)  # 1 hour

        # Mock grep output with LOOP entries (fast loops - 10 seconds each)
        # Format: "LOOP:  CPU time  real_time"
        grep_output = 'LOOP:  cpu time    4.0: real time    4.0\nLOOP:  cpu time    4.0: real time    4.0'
        mock_transport.exec_command_wait_return = (0, grep_output, '')
        mock_transport.get_attribute_return = MockFileStatResult(mtime=time.time())

        # Call function
        result = monitor_loop_time(mock_node, mock_transport)

        # Assertions - should be fine as 12.3s < 3600s/10 = 360s
        assert result is None

    def test_monitor_loop_time_slow_loops(self, mock_transport):
        """Test monitor_loop_time with slow electronic loops."""
        # Setup node with walltime
        mock_node = MockNode(walltime_limit=3600)  # 1 hour

        # Mock grep output with LOOP entries (slow loops - 400 seconds each)
        grep_output = 'LOOP:  cpu time    400.0: real time    400\nLOOP:  cpu time    400.0: real time    400.0'
        mock_transport.exec_command_wait_return = (0, grep_output, '')
        mock_transport.get_attribute_return = MockFileStatResult(mtime=time.time())

        # Call function (default minimum_electronic_loops=10)
        result = monitor_loop_time(mock_node, mock_transport)

        # Assertions - should detect slow loops as 450.3s > 3600s/10 = 360s
        assert result is not None
        assert 'Less than 10 electronic loop can be run' in result
        assert '400.00 seconds' in result

    def test_monitor_loop_time_custom_minimum_loops(self, mock_transport):
        """Test monitor_loop_time with custom minimum electronic loops."""
        # Setup node with walltime
        mock_node = MockNode(walltime_limit=3600)  # 1 hour

        # Mock grep output with LOOP entries (200 seconds each)
        grep_output = 'LOOP:  cpu time    400.0: real time    400\nLOOP:  cpu time    400.0: real time    400.0'
        mock_transport.exec_command_wait_return = (0, grep_output, '')
        mock_transport.get_attribute_return = MockFileStatResult(mtime=time.time())

        # Call function with minimum_electronic_loops=20
        result = monitor_loop_time(mock_node, mock_transport, minimum_electronic_loops=20)

        assert result is not None
        assert 'Less than 20 electronic loop can be run' in result
        assert '400.00 seconds' in result

    def test_monitor_loop_time_stalled_calculation(self, mock_transport, monkeypatch):
        """Test monitor_loop_time detecting stalled calculation."""
        # Setup time mock - current time is much later than file modification time
        current_time = 10000
        file_mtime = 8000  # File last modified 2000 seconds ago
        monkeypatch.setattr(time, 'time', lambda: current_time)

        # Setup node with walltime
        mock_node = MockNode(walltime_limit=7200)  # 2 hour

        # Mock grep output with LOOP entries (100 seconds each loop)
        grep_output = 'LOOP:  cpu time    100.0: real time    100\nLOOP:  cpu time    100.0: real time    100.0'
        mock_transport.exec_command_wait_return = (0, grep_output, '')
        mock_transport.get_attribute_return = MockFileStatResult(mtime=file_mtime)

        # Call function (default patience_num_loops=5, patience_minimum_time=1800)
        result = monitor_loop_time(mock_node, mock_transport)

        assert result is not None
        assert 'file is more than 500.00 seconds ago' in result

    def test_monitor_loop_time_recent_update(self, mock_transport, monkeypatch):
        """Test monitor_loop_time with recent file update."""
        # Setup time mock - file was modified recently
        current_time = 10000
        file_mtime = 9900  # File last modified 100 seconds ago
        monkeypatch.setattr(time, 'time', lambda: current_time)

        # Setup node with walltime
        mock_node = MockNode(walltime_limit=7200)  # 2 hour

        # Mock grep output with LOOP entries (50 seconds each loop)
        grep_output = 'LOOP:  cpu time    50.5: real time    50.5\nLOOP:  cpu time    50.5: real time    50.5'
        mock_transport.exec_command_wait_return = (0, grep_output, '')
        mock_transport.get_attribute_return = MockFileStatResult(mtime=file_mtime)

        # Call function
        result = monitor_loop_time(mock_node, mock_transport)

        # Assertions
        # elapsed = 10000 - 9900 = 100s
        # last_loop_time * patience_num_loops = 50.5 * 5 = 252.5s > 100s (not stalled)
        assert result is None

    def test_monitor_loop_time_minimum_patience_not_met(self, mock_transport, monkeypatch):
        """Test monitor_loop_time when minimum patience time is not met."""
        # Setup time mock
        current_time = 10000
        file_mtime = 8500  # File last modified 1500 seconds ago
        monkeypatch.setattr(time, 'time', lambda: current_time)

        # Setup node with walltime
        mock_node = MockNode(walltime_limit=7200)  # 2 hour

        # Mock grep output with LOOP entries
        grep_output = 'LOOP:  cpu time    100.5: real time    100.5\nLOOP:  cpu time    100.5: real time    100.5'
        mock_transport.exec_command_wait_return = (0, grep_output, '')
        mock_transport.get_attribute_return = MockFileStatResult(mtime=file_mtime)

        # Call function
        result = monitor_loop_time(mock_node, mock_transport)

        assert result is None
