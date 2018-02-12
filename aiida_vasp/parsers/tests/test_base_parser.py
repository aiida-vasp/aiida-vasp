"""Unit tests for BaseParser"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
import os

import pytest
from aiida.common.folders import SandboxFolder

from aiida_vasp.parsers.base import BaseParser
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC


@pytest.fixture()
def base_parser(vasp_nscf_and_ref):
    calc, _ = vasp_nscf_and_ref
    return BaseParser(calc)


@ONLY_ONE_CALC
def test_parse_with_retrieved(ref_retrieved_nscf, base_parser):
    assert base_parser.parse_with_retrieved({'retrieved': ref_retrieved_nscf})
    assert base_parser.out_folder == ref_retrieved_nscf


@ONLY_ONE_CALC
def test_result(base_parser):
    """Test result function returns success and new nodes"""
    base_parser.add_node('test', 'result')
    result, nodes = base_parser.result(True)
    assert result
    assert nodes == [('test', 'result')]

    result, nodes = base_parser.result(False)
    assert not result
    assert nodes == [('test', 'result')]


@ONLY_ONE_CALC
def test_add_node(base_parser):
    base_parser.add_node('test', 'result')
    assert base_parser.new_nodes['test'] == 'result'


@ONLY_ONE_CALC
def test_get_file(base_parser, ref_retrieved_nscf):
    """Test getting a retrieved output file."""
    base_parser.parse_with_retrieved({'retrieved': ref_retrieved_nscf})
    assert os.path.isfile(base_parser.get_file('OUTCAR'))
    assert os.path.exists(base_parser.get_file('OUTCAR'))
    assert base_parser.get_file('NonExistent') is None
