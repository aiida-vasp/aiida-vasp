"""Test the standard stream parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy as np

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.parsers.node_composer import NodeComposer


def test_stream_parser(fresh_aiida_env, stream_parser):
    """Test that the stream parser works."""
    import re  # pylint: disable=import-outside-toplevel
    errors = stream_parser.errors
    assert errors
    assert errors[0].kind == 'ERROR'
    assert errors[0].regex == re.compile('internal error in subroutine IBZKPT')
    assert len(stream_parser.warnings) == 0
