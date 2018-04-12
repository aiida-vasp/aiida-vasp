"""Module containing various decorators implementing delegate types."""
from functools import update_wrapper


def delegate():
    """
    Get a decorator adding attributes to add or remove functions to a list of functions.

    :param return_type: Decides how many items will be returned. Must be one of 'void',
                        'single' and 'list', meaning no, one or many return values.

    When the decorated function is called, all functions in the list will be called. It
    will return a list of all the return values of the subscribed functions. If the
    delegate does not have any subscribers the code from the delegate method itself will
    be executed. The basic usage is outlined in the corresponding utils.tests.test_delegates.
    For a fully fleshed out example look at the VaspParser and the BaseFileParser. The idea is
    that there is a class e.g.

    class VaspParser():

        @delegate()
        def parse_quantity(self, quantity, inputs):
            return None

    managing a number of instances of classes e.g.

    class FileParser():

        def __init__(self, cls):
            cls.parse_quantity.add_listener(self.parse_file)

        def parse_file(self, quantity, inputs):

            if can_parse(quantity):
               return parse(quantity, inputs)

            return None

    The FileParsers can be instantiated somewhere in the code by giving them an instance
    of the VaspParser i.e.

    vasp_parser = VaspParser()
    FileParser(vasp_parser)

    and they will subscribe their parse_file method to the parse_quantity delegate of the VaspParser.
    The VaspParser can then parse a given quantity by calling its parse_quantity method

    self.parse_quantity(inputs)

    which will then call all the methods subscribed to it. The advantage is that the
    VaspParser itself, does not need to keep track of all the FileParser objects. In addition
    it does not need to decide, which of the FileParsers to call. This logic can be moved to
    the FileParsers, who may decide whether they can parse that quantity based on the input.
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
            """Wrap the method to be delegated."""
            results = []
            for func in meth.listeners:
                results.append(func(*args[1:], **kwargs))
            for result in results:
                if result:
                    return result
            return meth(*args, **kwargs)

        update_wrapper(wrapper, meth)
        return wrapper

    return decorator
