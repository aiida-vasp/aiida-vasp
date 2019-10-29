"""
Extensions of dictionaries.

---------------------------
Extensions of Pythons standard dict as well as Aiida's AttributeDict.
"""

from aiida.common.extendeddicts import AttributeDict


class DictWithAttributes(AttributeDict):
    """
    Extension of the AttributeDict from Aiida.common.

    This class internally stores values in a dictionary, but exposes
    the keys also as attributes, i.e. asking for attrdict.key
    will return the value of attrdict['key'] and so on.

    If the key is not in the dict a default value will be returned.
    """

    def __getattr__(self, attr):
        """Read a key as an attribute. Return a Default value on missing key."""
        return self.get(attr)

    def __setattr__(self, attr, value):
        """Set a key as an attribute."""
        self[attr] = value
