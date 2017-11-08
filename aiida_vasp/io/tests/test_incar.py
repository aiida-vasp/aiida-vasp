"""Test the Incar io interface"""
import re
import pytest

from aiida_vasp.utils.fixtures.testdata import data_path, read_file
from aiida_vasp.io.incar import IncarIo


@pytest.fixture()
def incar_dict():
    incar_dict = {
        'encut': 350,
        'sigma': .5e-1,
        'lreal': False,
        'prec': 'Accurate'
    }
    return incar_dict


def test_read_incar():
    incar_path = data_path('phonondb', 'INCAR')
    incar_io = IncarIo(file_path=incar_path)
    assert incar_io.incar_dict['prec'] == 'Accurate'
    assert incar_io.incar_dict['ibrion'] == -1
    assert incar_io.incar_dict['encut'] == 359.7399
    assert incar_io.incar_dict['lreal'] == False


def test_from_dict(incar_dict):
    incar_io = IncarIo(incar_dict=incar_dict)
    assert str(
        incar_io
    ) == 'ENCUT = 350\nLREAL = .False.\nPREC = Accurate\nSIGMA = 0.05\n'


def test_write_incar(tmpdir, incar_dict):
    incar_io = IncarIo(incar_dict=incar_dict)
    tempfile = str(tmpdir.join('INCAR'))
    incar_io.store(tempfile)
    assert read_file(path=tempfile) == str(incar_io)
