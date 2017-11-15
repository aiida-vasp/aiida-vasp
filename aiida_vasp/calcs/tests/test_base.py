"""Unit tests for aiida_vasp.calcs.base"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import pytest
from aiida.common.exceptions import ValidationError

from aiida_vasp.utils.fixtures import aiida_env, fresh_aiida_env

TYPE = 'parameter'
DOC = 'input_parameters'
PARAM = 'kind'


def linkname(param):
    return 'linkname_{}'.format(param)


@pytest.mark.usefixtures('fresh_aiida_env')
def test_input_get_dict():
    """get_dict must have same format as documented for _use_methods item"""
    from aiida.orm import DataFactory
    from aiida_vasp.calcs.base import Input

    input_obj = Input(types=TYPE, doc=DOC, param=PARAM, ln=linkname)
    use_spec = input_obj.get_dict()
    assert use_spec['valid_types'] == (DataFactory(TYPE),)
    assert use_spec['additional_parameter'] == PARAM
    assert use_spec['linkname'] == linkname
    assert use_spec['docstring'] == DOC


@pytest.mark.usefixtures('fresh_aiida_env')
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


@pytest.fixture()
def base_calc():
    from aiida_vasp.calcs.base import VaspCalcBase
    return VaspCalcBase()


def test_check_inputs_fail(fresh_aiida_env, base_calc):
    with pytest.raises(ValidationError):
        base_calc.check_input(base_calc.get_inputs_dict(), 'code')
    with pytest.raises(ValidationError):
        base_calc.verify_inputs(base_calc.get_inputs_dict())


def test_check_inputs_success(fresh_aiida_env, base_calc):
    from aiida.orm import Code
    base_calc.use_code(Code())
    assert base_calc.check_input(base_calc.get_inputs_dict(), 'code')
    assert base_calc.verify_inputs(base_calc.get_inputs_dict())
