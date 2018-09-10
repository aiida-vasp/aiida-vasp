"""Inheritance tools."""
from textwrap import dedent


def update_docstring(method_name, content, append=True):
    r"""
    Update docstring of (an inherited) class method.

    For subclasses that use hooks to change behaviour of superclass methods.

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

    def wrapper(cls):
        """Update the method docstring and return the class."""
        if append:
            old_doc = getattr(cls, method_name).im_func.func_doc
            getattr(cls, method_name).im_func.func_doc = dedent(old_doc) + dedent(content)
        else:
            getattr(cls, method_name).im_func.func_doc = content
        return cls

    return wrapper
