"""Test the Incar io interface"""
# pylint: disable=redefined-outer-name
from collections import OrderedDict

import pytest

from aiida_vasp.utils.fixtures.testdata import data_path, read_file
from aiida_vasp.io.incar import IncarIo, IncarItem, IncarParamParser


@pytest.fixture()
def incar_dict():
    """Create a mapping of mixed case names to mixed parsed / unparsed values."""
    incar_dict = OrderedDict([('encut', 350), ('Sigma', '.5e-1 comment'), ('lreal', False), ('PREC', 'Accurate')])
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
    assert isinstance(incar_dict['bmix'], int)

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
