"""Test the Incar io interface"""
# pylint: disable=redefined-outer-name
from collections import OrderedDict

import pytest

from aiida_vasp.utils.fixtures.testdata import data_path, read_file
from aiida_vasp.io.incar import IncarIo, IncarItem, IncarParamParser, IncarParser
from aiida_vasp.utils.aiida_utils import get_data_class


@pytest.fixture()
def incar_dict():
    """Create a mapping of mixed case names to mixed parsed / unparsed values."""

    incar_dict = OrderedDict([('encut', 350), ('Sigma', '.5e-1 comment'), ('lreal', False), ('PREC', 'Accurate')])
    return incar_dict


@pytest.fixture()
def incar_dict_example():
    """Create a example dict."""

    incar_dict = {'encut': 350, 'Sigma': '.5e-1 #comment', 'lreal': False, 'PREC': 'Accurate'}
    return incar_dict


@pytest.mark.incar
def test_read_incar():
    """Read an INCAR file and test that some of the keys are read correctly."""

    incar_path = data_path('phonondb', 'INCAR')
    incar_io = IncarIo(file_path=incar_path)
    incar_dict = incar_io.get_dict()
    assert incar_dict['prec'] == 'Accurate'
    assert incar_dict['ibrion'] == -1
    assert incar_dict['encut'] == 359.7399
    assert incar_dict['lreal'] is False


@pytest.mark.incar
def test_example_incar():
    """Read a pathological case of an INCAR file (top level example from VASP docs)."""

    incar_path = data_path('incar_set', 'INCAR.copper_srf')
    incar_io = IncarIo(file_path=incar_path)
    incar_dict = incar_io.get_dict()
    assert incar_dict['system'] == 'Copper surface calculation'

    assert incar_dict['istart'] == 0
    assert isinstance(incar_dict['istart'], int)

    assert incar_dict['encut'] == 200.01
    assert isinstance(incar_dict['encut'], float)

    assert incar_dict['bmix'] == 2.0
    assert isinstance(incar_dict['bmix'], float)

    assert incar_dict['nelmin'] == 0
    assert incar_dict['nelmdl'] == 3


@pytest.mark.incar
def test_from_dict(incar_dict):
    incar_io = IncarIo(incar_dict=incar_dict)
    ref_str = '\n'.join(sorted(['ENCUT = 350', 'SIGMA = 0.05', 'LREAL = .False.', 'PREC = Accurate']))
    assert str(incar_io) == ref_str


@pytest.mark.incar
def test_from_string():
    """Test reading from string."""

    test_str = 'TRUE = .True\nFALSE=.f.'
    incar_io = IncarIo()
    incar_io.read_string(test_str)
    incar_dict = incar_io.get_dict()
    assert incar_dict.pop('true') is True
    assert incar_dict.pop('false') is False
    assert not incar_dict


@pytest.mark.incar
def test_write_incar(tmpdir, incar_dict):
    """Test writing and INCAR file from an IncarIo object."""

    incar_io = IncarIo(incar_dict=incar_dict)
    tempfile = str(tmpdir.join('INCAR'))
    incar_io.write(tempfile)
    assert read_file(path=tempfile) == str(incar_io)


@pytest.mark.incar
def test_incar_item():
    """Test the incar item class used to write to file"""

    test_str = 'ENCUT = 350 # test comment'
    item = IncarItem('encut', 350, '# test comment')
    assert item.name == 'ENCUT'
    assert item.value == 350
    assert item.comment == 'test comment'
    assert str(item) == test_str

    item = IncarItem(name='encut', value=350, comment='# test comment')
    assert item.name == 'ENCUT'
    assert item.value == 350
    assert item.comment == 'test comment'
    assert str(item) == test_str

    item = IncarItem.from_string(test_str)
    assert str(item) == test_str


@pytest.mark.incar
def test_parser():
    """Test the parser with a pathological string example."""

    test_string = '''TRUE = .True.
    FALSE=.False. this is a comment; FLOAT\t=\t1.45e-03
    INT = 45  # endline comment; may contain '#' and ';' NOPARAM = this is not a parameter
    LIST = 1 2 -33 5.6
    '''
    parsed = IncarParamParser.parse_string(test_string)
    assert parsed['true'] is True
    assert parsed['false'] is False
    assert parsed['float'] == 1.45e-3
    assert parsed['list'] == [1, 2, -33, 5.6]
    assert parsed['int'] == 45
    assert 'noparam' not in parsed


@pytest.mark.incar
def test_parser_read_parsevasp():
    """Test to read a INCAR file from parsevasp."""

    path = data_path('phonondb', 'INCAR')
    parser = IncarParser(file_path=path)
    result = parser.get_quantity('incar', {})
    assert isinstance(result['incar'], get_data_class('parameter'))
    incar = result['incar'].get_dict()
    assert incar['prec'] == 'Accurate'
    assert incar['ibrion'] == -1
    assert incar['encut'] == 359.7399
    assert incar['lreal'] is False


@pytest.mark.incar
def test_parser_read_doc_parsevasp():
    """
    Read example INCAR from VASP documentation.

    Using parsevasp.This test should fail as the comment line does
    not start with hashtag.

    """

    path = data_path('incar_set', 'INCAR.copper_srf')
    try:
        IncarParser(file_path=path)
    except AssertionError:
        pass


@pytest.mark.incar
def test_parser_dict_parsevasp(incar_dict_example):
    """
    Pass a dict to the INCAR parser.

    Using parsevasp. Should fail, since passing of dict in
    the interface is not implemented yet.

    """

    try:
        IncarParser(incar_dict=incar_dict_example)
    except AttributeError:
        pass


@pytest.mark.incar
def test_parser_string_parsevasp():
    """
    Pass a string to the INCAR parser.

    Using parsevasp. Should fail, since passing of string in
    the interface is not implemented yet.

    """

    test_string = 'TRUE = .True.\nFalse=.false.'
    try:
        IncarParser(incar_string=test_string)
    except AttributeError:
        pass


@pytest.mark.incar
def test_parser_write_parser(tmpdir, incar_dict_example):
    """Test writing an INCAR from a dict, read and compare."""

    # create ParameterData instances
    incar_params = get_data_class('parameter')(dict=incar_dict_example)
    assert isinstance(incar_params, get_data_class('parameter'))
    parser = IncarParser(data=incar_params)
    result = parser.get_quantity('incar', {})
    assert isinstance(result['incar'], get_data_class('parameter'))
    # now write
    temp_file = str(tmpdir.join('INCAR'))
    parser.write(temp_file)
    # read again
    parser_reparse = IncarParser(file_path=temp_file)
    result = parser_reparse.get_quantity('incar', {})
    assert isinstance(result['incar'], get_data_class('parameter'))
    comp_dict = {'encut': 350, 'sigma': 0.05, 'lreal': False, 'prec': 'Accurate'}
    assert str(sorted(result['incar'].get_dict())) == str(sorted(comp_dict))
