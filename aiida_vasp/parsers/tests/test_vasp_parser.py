"""Unittests for VaspParser"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
# pylint: disable=protected-access,unused-variable,too-few-public-methods

import pytest

from aiida_vasp.io.parser import BaseFileParser
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC


class TestFileParser(object):

    PARSABLE_ITEMS = {
        'quantity_with_alternatives': {
            'inputs': [],
            'parsers': ['DUMMY'],
            'nodeName': 'structure',
            'prerequisites': [],
        },
        'quantity1': {
            'inputs': [],
            'parsers': ['CONTCAR'],
            'nodeName': '',
            'is_alternative': 'quantity_with_alternatives',
            'prerequisites': []
        },
        'quantity2': {
            'inputs': [],
            'parsers': ['CONTCAR'],
            'nodeName': '',
            'prerequisites': ['quantity1']
        },
        'quantity3': {
            'inputs': [],
            'parsers': ['CONTCAR'],
            'nodeName': '',
            'is_alternative': 'non_existing_quantity',
            'prerequisites': ['quantity_with_alternatives']
        },
    }


@pytest.fixture
def vasp_parser_with_test(vasp_nscf_and_ref, ref_retrieved_nscf):
    """Fixture providing a VaspParser instance coupled to a VaspCalculation."""
    from aiida.orm.data.parameter import ParameterData
    vasp_calc, _ = vasp_nscf_and_ref
    vasp_calc.use_settings(ParameterData(dict={'parser_settings': {'add_quantity_with_alternatives': True, 'add_quantity2': True}}))
    parser = vasp_calc.get_parserclass()(vasp_calc)
    parser.add_file_parser('test_parser', {'parser_class': TestFileParser, 'is_critical': False})
    success, outputs = parser.parse_with_retrieved({'retrieved': ref_retrieved_nscf})
    return parser


@ONLY_ONE_CALC
def test_parsable_quantities(vasp_parser_with_test):
    """Check whether parsable quantities are set as intended."""
    parser = vasp_parser_with_test
    parsable_quantities = parser._parsable_quantities
    for quantity in TestFileParser.PARSABLE_ITEMS:
        assert quantity in parsable_quantities

    print parsable_quantities['quantity1']
    assert parsable_quantities['quantity1'].have_files
    assert parsable_quantities['quantity1'].is_parsable
    assert not parsable_quantities['quantity_with_alternatives'].have_files
    assert parsable_quantities['quantity2'].is_parsable
    assert not parsable_quantities['quantity3'].is_parsable
    assert 'non_existing_quantity' in parsable_quantities


@ONLY_ONE_CALC
def test_quantities_to_parse(vasp_parser_with_test):
    """Check if quantities are added to quantities to parse correctly."""
    parser = vasp_parser_with_test
    quantities_to_parse = parser._quantities_to_parse
    # after parse_with_retrieved quantities_to_parse will be empty, so we
    # have to set it one more time.
    parser._check_and_validate_settings()
    assert 'quantity2' in quantities_to_parse
    assert 'quantity_with_alternatives' not in quantities_to_parse
    assert 'quantity1' in quantities_to_parse
