"""Unit tests for mock-vasp command"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import subprocess
import six

def test_mock_vasp():
    """A simple test to verify that mock-vasp can run."""
    if six.PY2:
        output = subprocess.check_output(["mock-vasp", "--help"])
    else:
        output = subprocess.run(["mock-vasp", "--help"], stdout=subprocess.PIPE)
