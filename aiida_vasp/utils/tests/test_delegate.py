"""Test the delegate decorator"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.delegates import delegate


def test_delegate():
    """Test the functionality of the delegate() wrapper"""

    @delegate()
    def my_delegate(inputs):
        return None

    def method(inputs):
        return inputs

    inputs = 'This is a test input'

    assert not my_delegate.listeners
    my_delegate.add_listener(method)
    assert my_delegate(inputs) == inputs
    my_delegate.remove_listener(method)
    assert my_delegate.listeners == []
