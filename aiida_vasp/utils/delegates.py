"""
Delegate types.

---------------
Module containing decorators and classes implementing delegate types.
"""
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


class Delegate(list):
    """A callable list, that will call every function inside the list and delegate the arguments."""

    def __call__(self, *args, **kwargs):
        """Call all the functions subscribed to this list and return their result."""
        results = []
        for func in self:
            results.append(func(*args, **kwargs))
        for result in results:
            if result:
                return result
        return {args[0]: None}

    def clear(self):
        """Clear the listener list, for Python > 3.3 this will be a built-in method."""
        del self[:]
