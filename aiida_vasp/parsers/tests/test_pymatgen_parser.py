"""Unittests for the PymatgenParser"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
import os

import pytest

from aiida_vasp.parsers.pymatgen_vasp import PymatgenParser
from aiida_vasp.utils.fixtures import *


def data_path(*args):
    """path to a test data file"""
    path = os.path.realpath(
        os.path.join(__file__, '../../../test_data', *args))
    assert os.path.exists(path)
    assert os.path.isabs(path)
    return path


@pytest.fixture()
def parse_result(aiida_env):
    """
    Result of parsing a retrieved calculation (emulated)

    1. create a calculation with parser settings
    2. create a parser with the calculation
    3. populate a fake retrieved folder and pass it to the parser
    """
    from aiida.orm import CalculationFactory, DataFactory
    calc = CalculationFactory('vasp.vasp')()
    calc.use_settings(
        DataFactory('parameter')(dict={
            'pymatgen_parser': {
                'parse_potcar_file': False
            }
        }))
    parser = PymatgenParser(calc=calc)
    retrieved = DataFactory('folder')()
    retrieved.add_path(data_path('phonondb', 'vasprun.xml'), '')
    result = parser.parse_with_retrieved({'retrieved': retrieved})
    assert result[0]
    return result


def test_kpoints_result(parse_result):
    from aiida.orm import DataFactory
    _, nodes = parse_result
    assert isinstance([n for n in nodes if n[0] == 'kpoints'][0][1],
                      DataFactory('array.kpoints'))


def test_structure_result(parse_result):
    from aiida.orm import DataFactory
    _, nodes = parse_result
    assert isinstance([n for n in nodes if n[0] == 'structure'][0][1],
                      DataFactory('structure'))
