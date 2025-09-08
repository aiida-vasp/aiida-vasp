"""Unit tests for mock-vasp command."""

import subprocess


def test_mock_vasp():
    """A simple test to verify that mock-vasp can run."""
    subprocess.check_output(['mock-vasp', '--help'])
