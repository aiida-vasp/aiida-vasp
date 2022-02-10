"""Test the standard stream parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.fixtures import *


@pytest.mark.parametrize(['stream_parser'], [('stdout/out',)], indirect=True)
def test_stream_parser(stream_parser):
    """Test that the stream parser works."""
    import re  # pylint: disable=import-outside-toplevel
    errors = stream_parser.errors
    assert errors
    assert errors[0].kind == 'ERROR'
    assert errors[0].regex == re.compile('internal error in subroutine IBZKPT')
    assert len(stream_parser.warnings) == 0
