"""Abstract base class for INCAR tag helper classes."""
import abc


class AbstractIncarParam(object):
    """
    Abstract base class for INCAR tag helper classes.

    Interface:

        * method ``validate()``: This method should validate that the current internal state represents a valid value for the tag.
        * method ``clean(incar_dict)``: Scan the other tags in ``incar_dict`` for problems related to this tag.
        * property ``info``: A string, containing additional information about the meaning of the current value.
        * property `` value``: The (python-) value that should be written into the incar dictionary
        * property ``name``: the name of the tag in all upper case
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def validate(self):
        pass

    @abc.abstractmethod
    def clean(self, incar_dict):
        pass

    @abc.abstractproperty
    def info(self):
        pass

    @abc.abstractproperty
    def value(self):
        pass

    @abc.abstractproperty
    def name(self):
        pass

    def __str__(self):
        return '{name} = {value} ({info})'.format(name=self.name, value=self.value, info=self.info)

    @property
    def param(self):
        return {self.name.lower(): self.value}
