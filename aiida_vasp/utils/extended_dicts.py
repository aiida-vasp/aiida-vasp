"""Extensions of Pythons standard dict as well as Aiida's extendedDicts."""

from aiida.common.extendeddicts import AttributeDict


class DictWithAttributes(AttributeDict):
    """
    Extension of the AttributeDict from Aiida.common.

    This class internally stores values in a dictionary, but exposes
    the keys also as attributes, i.e. asking for attrdict.key
    will return the value of attrdict['key'] and so on.

    If the key is not in the dict a default value will be returned. Every
    key, that is not defined in _DEFAULT_VALUES, will simply return None.
    """

    def __init__(self, init=None, defaults=None):
        """
        Possibly set the initial values of the dictionary from an external dictionary init.

        Note that the attribute-calling syntax will work only 1 level deep.
        """
        super(DictWithAttributes, self).__init__(init)
        if defaults is None:
            defaults = {}
        self._default_values = defaults

    def __getattr__(self, attr):
        """Read a key as an attribute. Return a Default value on missing key."""
        return self.get(attr)

    def __setattr__(self, attr, value):
        """Set a key as an attribute."""
        self[attr] = value

    def get(self, attr):
        default = self['_default_values'].get(attr, None)
        if attr in self:
            return self[attr]

        return default
