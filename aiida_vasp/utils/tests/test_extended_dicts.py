"""Test extended dicts."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import,no-member

import pytest

from aiida.common.extendeddicts import AttributeDict
from aiida.plugins import DataFactory

from aiida_vasp.utils.extended_dicts import update_nested_dict
from aiida_vasp.utils.extended_dicts import delete_keys_from_dict


@pytest.fixture
def init_dict1():
    """Fixture for dictionary one."""
    dct = AttributeDict()
    dct.dct1 = AttributeDict()
    dct.dct1.test = 1.0
    dct.dct1.test2 = 'string1'
    dct.dct1.test3 = [1.0]
    dct.dct1.test4 = {'key1': 'value1'}

    return dct


@pytest.fixture
def init_dict2():
    """Fixture for dictionary two."""
    dct = AttributeDict()
    dct.dct1 = AttributeDict()
    dct.dct1.test = 2.0
    dct.dct1.test2 = 'string2'
    dct.dct1.test3 = [2.0]
    dct.dct1.test4 = {'key2': 'value2'}

    return dct


def test_nested_dict_update(init_dict1, init_dict2):
    """Test a nested dict update."""
    dict1 = AttributeDict(init_dict1)
    dict2 = AttributeDict(init_dict2)
    update_nested_dict(dict1, dict2)
    test = AttributeDict(
        {'dct1': AttributeDict({
            'test': 2.0,
            'test2': 'string2',
            'test3': [2.0],
            'test4': {
                'key1': 'value1',
                'key2': 'value2'
            }
        })})
    assert test == dict1
    dict1 = AttributeDict(init_dict1)
    del dict2.dct1.test
    update_nested_dict(dict1, dict2)
    test.dct1.test = 1.0
    assert test == dict1


def test_nested_dict_delete(init_dict1):
    """Test nested dict delete."""
    dict1 = AttributeDict(init_dict1)
    delete_keys_from_dict(dict1, 'dct1.test')
    test = AttributeDict({'dct1': AttributeDict({'test2': 'string1', 'test3': [1.0], 'test4': AttributeDict({'key1': 'value1'})})})
    assert test == dict1
    delete_keys_from_dict(dict1, 'dct1.test5')
    assert test == dict1
    delete_keys_from_dict(dict1, 'dct2')
    assert test == dict1
