"""Module containing decorators and classes implementing delegate types."""


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
