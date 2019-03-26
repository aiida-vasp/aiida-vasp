
"""Test aiida_utils functionss."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import,no-member

import pytest
from importlib import import_module
from aiida_vasp.utils.aiida_utils import get_data_node, get_data_class, dbenv, BASIC_DATA_TYPES, get_current_user
from aiida_vasp.utils.fixtures.environment import aiida_env


def test_get_current_user(aiida_env):
    """Assert that get_current_user returns a user in all tested aiida versions."""
    user = get_current_user()
    assert user.first_name
    assert user.last_name
    assert user.email


def test_get_data_class(aiida_env):
    """Make sure the get_data_class accept valid types."""
    for data_type in BASIC_DATA_TYPES:
        data_type_class = get_data_class(data_type)
        the_module_ref =  __import__('aiida.orm', fromlist=[data_type.capitalize()])
        aiida_data_type_class = getattr(the_module_ref, data_type.capitalize())
        assert data_type_class == aiida_data_type_class

    with pytest.raises(KeyError) as e_info:
        get_data_class('garbage')

def test_get_data_node(aiida_env):
    """Make sure the get_data_node returns objects for the basic data types."""
    for data_type in BASIC_DATA_TYPES:
        the_module_ref =  __import__('aiida.orm', fromlist=[data_type.capitalize()])
        aiida_data_type_class = getattr(the_module_ref, data_type.capitalize())
        if data_type == 'bool':
            data_node = get_data_node(data_type, True)
            aiida_data_node = aiida_data_type_class(True)
            assert data_node.value == aiida_data_node.value
        if data_type == 'int':
            data_node = get_data_node(data_type, 1)
            aiida_data_node = aiida_data_type_class(1)
            assert data_node.value == aiida_data_node.value
        if data_type == 'float':
            data_node = get_data_node(data_type, 1.0)
            aiida_data_node = aiida_data_type_class(1.0)
            assert data_node.value == aiida_data_node.value
        if data_type == 'str':
            data_node = get_data_node(data_type, '')
            aiida_data_node = aiida_data_type_class('')
            assert data_node.value == aiida_data_node.value
        if data_type == 'list':
            data_node = get_data_node(data_type, list=[])
            aiida_data_node = aiida_data_type_class(list=[])
            assert data_node.get_list() == aiida_data_node.get_list()
        if data_type == 'dict':
            data_node = get_data_node(data_type, dict={})
            aiida_data_node = aiida_data_type_class(dict={})
            assert set(data_node.get_dict()) == set(aiida_data_node.get_dict())

    with pytest.raises(KeyError) as e_info:
        get_data_node('garbage', True)
    
