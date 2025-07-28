"""
Delegate types.

Module containing decorators and classes implementing delegate types.
"""

from __future__ import annotations

from functools import update_wrapper
from typing import Any, Callable, TypeVar

F = TypeVar('F', bound=Callable[..., Any])


def delegate_method_kwargs(prefix: str = '_init_with_') -> Callable[[F], F]:
    """
    Get a kwargs delegating decorator.

    :params prefix: (str) common prefix of delegate functions
    """

    def decorator(meth: F) -> F:
        """Decorate a class method to delegate kwargs."""

        def wrapper(*args: Any, **kwargs: Any) -> Any:
            for kwarg, value in kwargs.items():
                getattr(args[0], prefix + kwarg)(value)
            meth(*args, **kwargs)

        update_wrapper(wrapper, meth)
        return wrapper

    return decorator
