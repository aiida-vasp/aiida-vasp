#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test fixtures"""

import json
import operator

import pytest


@pytest.fixture
# pylint: disable=protected-access
def test_name(request):
    """Returns module_name.function_name for a given test"""
    return request.module.__name__ + '/' + request._parent_request._pyfuncitem.name


@pytest.fixture
def compare_data(request, test_name_):
    """
    Returns a function which either saves some data to a file or (if that file exists already)

    compares it to pre-existing data using a given comparison function.
    """

    def inner(compare_fct, data, tag=None):
        full_name = test_name_ + (tag or '')
        val = request.config.cache.get(full_name, None)
        if val is None:
            request.config.cache.set(full_name, json.loads(json.dumps(data)))
            raise ValueError('Reference data does not exist.')
        else:
            assert compare_fct(val, json.loads(
                json.dumps(data)))  # get rid of json-specific quirks

    return inner


@pytest.fixture
def compare_equal(compare_data_):
    return lambda data, tag=None: compare_data_(operator.eq, data, tag)
