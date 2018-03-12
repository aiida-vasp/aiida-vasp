"""Utilities for working with aiida in general"""
from functools import wraps


def load_dbenv_if_not_loaded(**kwargs):
    """Load dbenv if necessary, run spinner meanwhile to show command hasn't crashed."""
    from aiida.backends.utils import load_dbenv, is_dbenv_loaded
    if not is_dbenv_loaded():
        load_dbenv(**kwargs)


def dbenv(function):
    """A function decorator that loads the dbenv if necessary before running the function."""

    @wraps(function)
    def decorated_function(*args, **kwargs):
        """Load dbenv if not yet loaded, then run the original function."""
        load_dbenv_if_not_loaded()
        return function(*args, **kwargs)

    return decorated_function


@dbenv
def get_data_node(data_type, *args, **kwargs):
    from aiida.orm import DataFactory
    return DataFactory(data_type)(*args, **kwargs)


@dbenv
def get_data_class(data_type):
    from aiida.orm import DataFactory
    return DataFactory(data_type)
