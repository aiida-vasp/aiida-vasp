"""
Inheritance tools.

This sets up inheritance of the docstrings for use when we inherit
classes. This makes it simple to add further details to or replace the
docstring present on the base class.
"""

from __future__ import annotations

from typing import Callable, TypeVar

ClassType = TypeVar('ClassType')


def update_docstring(
    method_name: str, content: str, append: bool = True
) -> Callable[[type[ClassType]], type[ClassType]]:
    r"""
    Update docstring of (an inherited) class method.

    For subclasses that use hooks to change behavior of superclass methods.

    :param method_name: Name of the method to update
    :param content: Content to add or replace
    :param append: If true, append to the docstring, else overwrite it entirely.

    Example::

        class Base(object):
            def method(self, **kwargs):
                '''Print the base_arg kwarg.'''
                print kwargs.pop('base_arg')
                self.process_additional(**kwargs)

            def process_additional(self, **kwargs):
                pass


        @update_docstring('method', '\n\nAlso print the sub_arg kwarg.', append=True)
        class Subclass(Base):
            def process_additional(self, **kwargs):
                print kwargs['sub_arg']


        @update_docstring('method', 'Print all kwargs.', append=False)
        class Allprinter(Base):
            def process_additional(self, **kwargs):
                for _, value in kwargs.items():
                    print value
    """

    def wrapper(cls: type[ClassType]) -> type[ClassType]:
        """Update the method docstring and return the class."""
        if append:
            getattr(cls, method_name).__func__.__doc__ += content
        else:
            getattr(cls, method_name).__func__.__doc__ = content
        return cls

    return wrapper
