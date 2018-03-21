"""Module containing various decorators implementing delegate types."""
from functools import update_wrapper


def delegate_method_kwargs(prefix='_init_with_'):
    """
    Get a kwargs delegating decorator.

    :params prefix: (str) common prefix of delegate functions
    """

    def decorator(meth):
        """Decorate a class method to delegate kwargs."""

        def wrapper(*args, **kwargs):
            for kwarg, value in kwargs.items():
                getattr(args[0], prefix + kwarg)(value)
            meth(*args, **kwargs)

        update_wrapper(wrapper, meth)
        return wrapper

    return decorator


def delegate():
    """
    Get a decorator adding attributes to add or remove functions to a list of functions.

    When the decorated function is called, all functions in the list will be called.
    """

    def decorator(meth):
        """Decorate a class method to delegate kwargs."""

        meth.listeners = []

        def add_listener(func):
            meth.listeners.append(func)

        def remove_listener(func):
            if func in meth.listeners:
                meth.listeners.remove(func)

        setattr(meth, 'add_listener', add_listener)
        setattr(meth, 'remove_listener', remove_listener)

        def wrapper(*args, **kwargs):
            for func in meth.listeners:
                func(*args, **kwargs)
            meth(*args, **kwargs)

        update_wrapper(wrapper, meth)
        return wrapper

    return decorator
