"""Unit tests for aiida_vasp.calcs.base"""
# pylint: disable=redefined-outer-name,unused-argument
import pytest
from aiida.utils.fixtures import FixtureManager

TYPE = 'parameter'
DOC = 'input_parameters'
PARAM = 'kind'


def linkname(param):
    return 'linkname_{}'.format(param)


@pytest.fixture
def load_dbenv():
    from aiida import load_dbenv
    load_dbenv()


@pytest.fixture
def create_aiida_env():
    manager = FixtureManager()
    manager.create_aiida_db()
    manager.create_profile()


@pytest.fixture
def input_obj(create_aiida_env):
    from aiida_vasp.calcs.base import Input

    return Input(types=TYPE, doc=DOC, param=PARAM, ln=linkname)


def test_input_get_dict(input_obj):
    """get_dict must have same format as documented for _use_methods item"""
    from aiida.orm import DataFactory
    use_spec = input_obj.get_dict()
    assert use_spec['valid_types'] == (DataFactory(TYPE), )
    assert use_spec['additional_parameter'] == PARAM
    assert use_spec['linkname'] == linkname
    assert use_spec['docstring'] == DOC


def test_input_type_seq():
    """make sure sequence of type strings is correctly converted"""
    from aiida.orm import DataFactory
    from aiida_vasp.calcs.base import Input
    type_1 = 'parameter'
    type_2 = 'structure'
    input_obj = Input(types=[type_1, type_2])
    type_1_cls = DataFactory(type_1)
    type_2_cls = DataFactory(type_2)
    assert input_obj.types == (type_1_cls, type_2_cls)
