"""Test the delegate decorator"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.delegates import delegate


def test_delegate():
    """Test the functionality of the delegate() wrapper"""

    @delegate()
    def my_delegate(inputs):
        pass

    def method(inputs):
        return inputs

    inputs = 'This is a test input'

    assert my_delegate.listeners is False
    my_delegate.add_listener(method)
    assert my_delegate(inputs)[0] == inputs
    my_delegate.remove_listener(method)
    assert my_delegate.listeners == []
