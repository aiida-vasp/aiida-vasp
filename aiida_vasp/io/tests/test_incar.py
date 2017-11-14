"""Test the Incar io interface"""
import re
import pytest
import py
from collections import OrderedDict

from aiida_vasp.utils.fixtures.testdata import data_path, read_file
from aiida_vasp.io.incar import IncarIo, IncarItem, IncarParamParser


@pytest.fixture()
def incar_dict():
    incar_dict = OrderedDict([('encut', 350), ('sigma', .5e-1),
                              ('lreal', False), ('prec', 'Accurate')])
    return incar_dict


def test_read_incar():
    incar_path = data_path('phonondb', 'INCAR')
    incar_io = IncarIo(file_path=incar_path)
    print incar_io.get_dict()
    assert incar_io.get_dict()['prec'] == 'Accurate'
    assert incar_io.get_dict()['ibrion'] == -1
    assert incar_io.get_dict()['encut'] == 359.7399
    assert incar_io.get_dict()['lreal'] == False


def test_from_dict(incar_dict):
    incar_io = IncarIo(incar_dict=incar_dict)
    assert str(
        incar_io
    ) == 'ENCUT = 350\nSIGMA = 0.05\nLREAL = .False.\nPREC = Accurate'


def test_write_incar(tmpdir, incar_dict):
    incar_io = IncarIo(incar_dict=incar_dict)
    tempfile = str(tmpdir.join('INCAR'))
    incar_io.store(tempfile)
    assert read_file(path=tempfile) == str(incar_io)


def test_incar_item():
    item = IncarItem('encut', 350, '# test comment')
    assert item.name == 'ENCUT'
    assert item.value == 350
    assert item.comment == 'test comment'
    assert str(item) == 'ENCUT = 350 # test comment'

    item = IncarItem(name='encut', value=350, comment='# test comment')
    assert item.name == 'ENCUT'
    assert item.value == 350
    assert item.comment == 'test comment'
    assert str(item) == 'ENCUT = 350 # test comment'


@pytest.mark.wip
def test_parser():
    test_string = '''TRUE = .True.
    FALSE=.False. this is a comment; FLOAT\t=\t1.45e-03
    LIST = 1 2 -33 5.6
    '''
    parsed = IncarParamParser.parse_string(test_string)
    print parsed
    assert parsed['true'] == True
    assert parsed['false'] == False
    assert parsed['float'] == 1.45e-3
    assert parsed['list'] == [1, 2, -33, 5.6]
