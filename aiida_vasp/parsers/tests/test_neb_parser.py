"""Unittests for VtstNebVaspParser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
# pylint: disable=protected-access,unused-variable,too-few-public-methods, import-outside-toplevel

import os

import numpy as np
import pytest

from aiida.common.links import LinkType
from aiida.plugins import CalculationFactory, ParserFactory

from aiida_vasp.parsers.content_parsers.base import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC, calc_with_retrieved, neb_calc_with_retrieved
from aiida_vasp.utils.fixtures.testdata import data_path


def _get_neb_vasp_parser(neb_calc_with_retrieved):
    """Return vasp parser before parsing"""
    settings_dict = {
        # 'ADDITIONAL_RETRIEVE_LIST': CalculationFactory('vasp.vasp')._ALWAYS_RETRIEVE_LIST,
        'parser_settings': {
            'add_image_forces': True
        }
    }
    file_path = str(os.path.abspath(os.path.dirname(__file__)) + '/../../test_data/neb')
    node = neb_calc_with_retrieved(file_path, settings_dict, 3)
    parser = ParserFactory('vasp.neb')(node)
    return parser, file_path, node


@pytest.fixture
def neb_parser_with_test(neb_calc_with_retrieved):
    """Fixture providing a VaspParser instance coupled to a VaspCalculation."""
    parser, file_path, node = _get_neb_vasp_parser(neb_calc_with_retrieved)
    success = parser.parse(retrieved_temporary_folder=file_path)
    try:
        yield parser
    finally:
        parser = ParserFactory('vasp.vasp')(node)


def test_neb_parser(neb_parser_with_test):
    """
    Test the neb parser
    """
    parser = neb_parser_with_test
    assert 'neb_misc' in neb_parser_with_test.outputs
    neb_misc = parser.outputs.neb_misc.get_dict()

    assert neb_misc['neb_data']['01']['neb_converged']
    assert neb_misc['neb_data']['01']['free_energy'] == -19.49164066

    assert neb_misc['neb_data']['03']

    # Check that notifications exists
    assert 'notifications' in parser.outputs['misc.image_01'].get_dict()
    assert 'notifications' in parser.outputs['misc.image_03'].get_dict()

    # Make sure structures are parsed as well
    assert 'structure.image_01' in parser.outputs

    # Check that the forces is parserd
    forces = parser.outputs.image_forces.get_array('forces_01')
    assert forces[0].tolist() == [0.008815, 0.005492, -0.000661]
    assert forces.shape == (4, 3)
