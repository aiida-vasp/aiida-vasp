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


@dbenv
def backend_obj_users():
    """Test if aiida accesses users through backend object."""
    backend_obj_flag = False
    try:
        from aiida.backends.utils import get_automatic_user  # pylint: disable=unused-variable,no-name-in-module
    except ImportError:
        backend_obj_flag = True
    return backend_obj_flag


@dbenv
def get_current_user():
    """Get current user backwards compatibly with aiida-core <= 0.12.1."""
    current_user = None
    if backend_obj_users():
        from aiida.orm.backend import construct_backend  # pylint: disable=no-name-in-module
        backend = construct_backend()
        current_user = backend.users.get_automatic_user()
    else:
        from aiida.backends.utils import get_automatic_user  # pylint: disable=no-name-in-module
        current_user = get_automatic_user()
    return current_user


def builder_interface(calc_cls):
    """Return the JobProcess or the JobCalculation class, depending on aiida version."""
    if hasattr(calc_cls, 'get_builder'):
        return True
    return False
