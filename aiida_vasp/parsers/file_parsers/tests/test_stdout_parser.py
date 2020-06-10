"""Test the STDOUT parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy as np

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.parsers.node_composer import NodeComposer


def test_stdout_error(fresh_aiida_env, stdout_parser):
    """Test parsing of stdout error"""

    stdout_parser.settings.nodes.update({'misc': {'type': 'dict', 'quantities': ['stdout_error'], 'link_name': 'my_custom_node'}})

    composer = NodeComposer(file_parsers=[stdout_parser])
    data_obj = composer.compose('dict', quantities=['stdout_error'])
    ref_class = get_data_class('dict')
    assert isinstance(data_obj, ref_class)
    data_dict = data_obj.get_dict()
    # test error
    test_dict = {'stdout_error': [{'shortname': 'zpotrf', 'message': 'Error in Lapack ZPOTRF', 'critical': True}]}

    assert data_dict == test_dict
