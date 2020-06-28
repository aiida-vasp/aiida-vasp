"""Tests for error handling module"""
import re
import io

import pytest
from aiida_vasp.parsers.error_handling import VaspError, VaspGeneralError, ErrorRecord, ErrorScanner

# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import


@pytest.fixture
def vasp_error_example():
    """Fixture for VaspError"""
    error = VaspError('BAD NEWS', shortname='bad', message='bad news error', critical=False, suggestion=None)
    return error


@pytest.fixture
def vasp_general_error_example():
    """Fixture for a GeneralError"""
    return VaspGeneralError()


def test_error_record():
    record = ErrorRecord(error='FAILURE', message='error', shortname='test_error', critical=False, suggestion=None)
    assert record.message == 'error'


def test_vasp_error(vasp_error_example):
    """Test the VaspError Class"""
    error = vasp_error_example
    assert isinstance(error.regex, re.Pattern)

    error.__repr__()
    erec = error.check_line('BAD NEWS: There is an error')
    assert erec is not None
    assert erec.shortname == error.shortname
    assert erec.critical is False
    assert erec.error is error


def test_vasp_general_error(vasp_general_error_example):
    """Test the VaspGeneralError Class"""
    error = vasp_general_error_example
    erec = error.check_line('FOO FAILURE:BAD DAY')
    assert erec is not None

    assert erec.critical is False
    assert erec.message == 'FAILURE:BAD DAY'
    assert erec.shortname == 'general'


@pytest.fixture
def error_stdout():
    """Fixture for a stream with errors"""

    content = """
BAD NEWS internal error in subroutine IBZKPT
FOO BAD NEWS VASP ERRORED
"""
    return io.StringIO(content)


def test_error_scanner(error_stdout):
    """Test ErrorScanner"""
    scanner = ErrorScanner(error_stdout)
    errors = list(scanner.get_errors())
    errors.sort(key=lambda x: x.shortname)
    assert errors[0].shortname == 'general'
    assert errors[0].message == 'BAD NEWS VASP ERRORED'
    assert errors[1].shortname == 'ibzkpt'
