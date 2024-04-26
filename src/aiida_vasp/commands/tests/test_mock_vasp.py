"""Unit tests for mock-vasp command."""

# pylint: disable=unused-import,unused-argument,redefined-outer-name
import subprocess


def test_mock_vasp():
    """A simple test to verify that mock-vasp can run."""
    subprocess.check_output(['mock-vasp', '--help'])
