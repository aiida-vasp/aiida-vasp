"""
Extensions of dictionaries.

---------------------------
Extensions of Pythons standard dict as well as Aiida's AttributeDict.
"""
import collections
from copy import deepcopy

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


def delete_keys_from_dict(dictionary, keys):
    """
    Delete a key from a nested dictionary.

    Extended to support somekey.someotherkey in case we need some restrictions on the nesting.
    """
    if not isinstance(keys, list):
        keylist = [keys]
    else:
        keylist = keys
    for key in keylist:
        nested_keys = key.strip().split('.')
        delete_nested_key(dictionary, nested_keys)


def delete_nested_key(dictionary, keys):
    """Delete the dictionary entry corresponding to a nested hierarchy of keys."""
    from collections.abc import MutableMapping  # pylint: disable=import-outside-toplevel
    from contextlib import suppress  # pylint: disable=import-outside-toplevel
    if keys and dictionary:
        element = keys[0]
        if element:
            value = dictionary.get(element)
            if len(keys) == 1:
                with suppress(KeyError):
                    del dictionary[element]
            else:
                if isinstance(value, MutableMapping):
                    delete_nested_key(value, keys[1:])


def update_nested_dict(dict1, dict2):
    """Updated a nested dictionary, where dict1 is updated with values in dict2."""
    for key, value in dict2.items():
        dict1_value = dict1.get(key)
        if isinstance(value, collections.Mapping) and isinstance(dict1_value, collections.Mapping):
            update_nested_dict(dict1_value, value)
        else:
            dict1[key] = deepcopy(value)


def find_key_in_dicts(dictionary, supplied_key):
    """Find a key in a nested dictionary."""
    for key, value in dictionary.items():
        if key == supplied_key:
            yield value
        elif isinstance(value, dict):
            for result in find_key_in_dicts(value, supplied_key):
                yield result
