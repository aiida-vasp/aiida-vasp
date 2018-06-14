"""Unittests for VaspParser"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
# pylint: disable=protected-access,unused-variable,too-few-public-methods

import pytest

from aiida_vasp.io.parser import BaseFileParser
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC


class ExampleFileParser(BaseFileParser):
    """Example FileParser class for testing VaspParsers functionality."""

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
            'nodeName': 'structure',
            'is_alternative': 'quantity_with_alternatives',
            'prerequisites': []
        },
        'quantity2': {
            'inputs': [],
            'parsers': ['_scheduler-stdout.txt'],
            'nodeName': 'trajectory',
            'is_alternative': 'trajectory',
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

    def __init__(self, *args, **kwargs):
        super(ExampleFileParser, self).__init__(*args, **kwargs)
        self._parsable_items = ExampleFileParser.PARSABLE_ITEMS
        self._parsable_data = {}

    def _parse_file(self, inputs):
        from aiida.orm.data.parameter import ParameterData
        result = {}
        for quantity in ExampleFileParser.PARSABLE_ITEMS:
            result[quantity] = ParameterData(dict={})
        return result


class ExampleFileParser2(BaseFileParser):
    """Example class for testing non unique quantity identifiers."""

    PARSABLE_ITEMS = {
        'quantity1': {
            'inputs': [],
            'parsers': ['CONTCAR'],
            'nodeName': '',
            'is_alternative': 'quantity_with_alternatives',
            'prerequisites': []
        },
    }

    def __init__(self, *args, **kwargs):
        super(ExampleFileParser2, self).__init__(*args, **kwargs)
        self._parsable_items = ExampleFileParser2.PARSABLE_ITEMS
        self._parsable_data = {}

    def _parse_file(self, inputs):
        from aiida.orm.data.parameter import ParameterData
        result = {}
        for quantity in ExampleFileParser.PARSABLE_ITEMS:
            result[quantity] = ParameterData(dict={})
        return result


@pytest.fixture
def vasp_parser_with_test(vasp_nscf_and_ref, ref_retrieved_nscf):
    """Fixture providing a VaspParser instance coupled to a VaspCalculation."""
    from aiida.orm.data.parameter import ParameterData
    vasp_calc, _ = vasp_nscf_and_ref
    vasp_calc.use_settings(ParameterData(dict={'parser_settings': {'add_quantity_with_alternatives': True, 'add_quantity2': True}}))
    parser = vasp_calc.get_parserclass()(vasp_calc)
    parser.add_file_parser('_scheduler-stdout.txt', {'parser_class': ExampleFileParser, 'is_critical': False})
    success, outputs = parser.parse_with_retrieved({'retrieved': ref_retrieved_nscf})
    try:
        yield parser
    finally:
        parser = vasp_calc.get_parserclass()(vasp_calc)


@ONLY_ONE_CALC
def test_quantities_to_parse(vasp_nscf_and_ref, ref_retrieved_nscf):
    """Check if quantities are added to quantities to parse correctly."""
    from aiida.orm.data.parameter import ParameterData
    vasp_calc, _ = vasp_nscf_and_ref
    vasp_calc.use_settings(ParameterData(dict={'parser_settings': {'add_quantity_with_alternatives': True, 'add_quantity2': True}}))
    parser = vasp_calc.get_parserclass()(vasp_calc)
    parser.add_file_parser('_scheduler-stdout.txt', {'parser_class': ExampleFileParser, 'is_critical': False})
    parser.out_folder = parser.get_folder({'retrieved': ref_retrieved_nscf})
    parser._set_parsable_quantities()
    parser._check_and_validate_settings()
    assert 'quantity2' in parser._quantities_to_parse
    assert 'quantity_with_alternatives' not in parser._quantities_to_parse
    assert 'quantity1' in parser._quantities_to_parse


@ONLY_ONE_CALC
def test_parsable_quantities(vasp_parser_with_test):
    """Check whether parsable quantities are set as intended."""
    parser = vasp_parser_with_test
    parsable_quantities = parser._parsable_quantities
    # Check whether all quantities from the added ExampleFileParser have been added.
    for quantity in ExampleFileParser.PARSABLE_ITEMS:
        assert quantity in parsable_quantities
    # Check whether quantities have been set up correctly.
    assert parsable_quantities['quantity1'].has_files
    assert parsable_quantities['quantity1'].is_parsable
    assert not parsable_quantities['quantity_with_alternatives'].has_files
    assert parsable_quantities['quantity2'].is_parsable
    assert not parsable_quantities['quantity3'].is_parsable
    # check whether the additional non existing quantity has been added. This is for cases,
    # where a quantity is an alternative to another main quantity, which has not been loaded.
    assert 'non_existing_quantity' in parsable_quantities


@ONLY_ONE_CALC
def test_node_link_names(vasp_parser_with_test):
    """Check whether an alternative quantity representing a node will be added with the correct linkname."""
    parser = vasp_parser_with_test

    print parser._parsers['_scheduler-stdout.txt'].parser._parsed_data
    assert 'quantity2' in parser._output_nodes
    print parser._output_nodes['quantity2']
    # 'quantity2' is alternative to 'trajectory', which is not going to be parsed here.
    assert 'output_trajectory' in parser._new_nodes


@ONLY_ONE_CALC
def test_quantity_uniqeness(vasp_parser_with_test):
    """Make sure non-unique quantity identifiers are detected."""
    parser = vasp_parser_with_test
    # Add a second ExampleFileParser that defines a quantity with the same identifier as the first one.
    parser.add_file_parser('another_test_parser', {'parser_class': ExampleFileParser2, 'is_critical': False})
    with pytest.raises(RuntimeError) as excinfo:
        parser._set_parsable_quantities()
    assert 'quantity1' in str(excinfo.value)
