"""Test the INCAR parser."""
# pylint: disable=redefined-outer-name, unused-wildcard-import, unused-argument, wildcard-import

import pytest
from pytest import raises
from aiida.common import InputValidationError

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.parsers.file_parsers.incar import IncarParser
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node


@pytest.fixture()
def incar_dict_example():
    """Create a example dict."""

    incar_dict = {'encut': 350, 'Sigma': '.5e-1 #comment', 'lreal': False, 'PREC': 'Accurate'}
    return incar_dict


@pytest.mark.incar
def test_parser_read(fresh_aiida_env):
    """Test to read a INCAR file."""

    path = data_path('phonondb', 'INCAR')
    parser = IncarParser(file_path=path)
    incar = parser.incar
    assert incar['prec'] == 'Accurate'
    assert incar['ibrion'] == -1
    assert incar['encut'] == 359.7399
    assert incar['lreal'] is False


@pytest.mark.incar
def test_parser_read_doc(fresh_aiida_env):
    """
    Read example INCAR from VASP documentation.

    Using parsevasp. Returned content should be none
    since parsevasp refuse to parse an INCAR where the
    comments does not start with hashtags.

    """

    path = data_path('incar', 'INCAR.copper_srf')
    parser = IncarParser(file_path=path)
    result = parser.incar
    assert result is None


@pytest.mark.incar
def test_parser_dict(fresh_aiida_env, incar_dict_example):
    """
    Pass a dict to the INCAR parser.

    Should return an AiiDA datastructure for dict.

    """

    parser = IncarParser(data=get_data_node('dict', dict=incar_dict_example))
    assert isinstance(parser.incar, get_data_class('dict'))


@pytest.mark.incar
def test_parser_string():
    """
    Pass a string to the INCAR parser.

    Should fail, since passing of string in
    the interface is not implemented yet.

    """

    test_string = 'LOPTICS = .True.\nAddgrid=.false.'
    with raises(AttributeError):
        IncarParser(incar_string=test_string)


@pytest.mark.incar
def test_write_parser(fresh_aiida_env, tmpdir, incar_dict_example):
    """Test writing an INCAR from a dict, read and compare."""

    # create AiiDA dictionary instance
    incar_params = get_data_class('dict')(dict=incar_dict_example)
    assert isinstance(incar_params, get_data_class('dict'))
    parser = IncarParser(data=incar_params)

    # now write
    temp_file = str(tmpdir.join('INCAR'))
    parser.write(temp_file)
    # read again
    parser_reparse = IncarParser(file_path=temp_file)
    result = parser_reparse.incar
    # compare
    comp_dict = {'encut': 350, 'sigma': 0.05, 'lreal': False, 'prec': 'Accurate'}
    assert str(sorted(result)) == str(sorted(comp_dict))

    # Test validation
    with_invalid = dict(incar_dict_example)
    with_invalid.update(foo='bar')
    incar_params = get_data_class('dict')(dict=with_invalid)
    parser = IncarParser(data=incar_params)
    with pytest.raises(InputValidationError):
        parser.write(temp_file)
