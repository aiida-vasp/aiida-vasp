"""
Extensions of dictionaries.

Extensions of Pythons standard dict as well as Aiida's AttributeDict.
"""

from __future__ import annotations

import collections.abc
from collections.abc import MutableMapping  # pylint: disable=import-outside-toplevel
from contextlib import suppress  # pylint: disable=import-outside-toplevel
from copy import deepcopy
from typing import Any

from aiida import orm
from aiida.common.extendeddicts import AttributeDict


class DictWithAttributes(AttributeDict):
    """
    Extension of the AttributeDict from Aiida.common.

    This class internally stores values in a dictionary, but exposes
    the keys also as attributes, i.e. asking for attrdict.key
    will return the value of attrdict['key'] and so on.

    If the key is not in the dict a default value will be returned.
    """

    def __getattr__(self, attr: str) -> Any:
        """Read a key as an attribute. Return a Default value on missing key."""
        return self.get(attr)

    def __setattr__(self, attr: str, value: Any) -> None:
        """Set a key as an attribute."""
        self[attr] = value


def delete_keys_from_dict(dictionary: dict[str, Any], keys: str | list[str]) -> None:
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


def delete_nested_key(dictionary: dict[str, Any], keys: list[str]) -> None:
    """Delete the dictionary entry corresponding to a nested hierarchy of keys."""

    if keys and dictionary:
        element = keys[0]
        if element:
            value = dictionary.get(element)
            if len(keys) == 1:
                with suppress(KeyError):
                    del dictionary[element]
            elif isinstance(value, MutableMapping):
                delete_nested_key(value, keys[1:])


def update_nested_dict(dict1: dict[str, Any], dict2: dict[str, Any], extend_list: bool = False) -> dict[str, Any]:
    """Updated a nested dictionary, where dict1 is updated with values in dict2."""
    for key, value in dict2.items():
        dict1_value = dict1.get(key)
        if isinstance(value, collections.abc.Mapping) and isinstance(dict1_value, collections.abc.Mapping):
            update_nested_dict(dict1_value, value)
        elif isinstance(value, list) and isinstance(dict1_value, list) and extend_list:
            dict1_value.extend(value)
        else:
            dict1[key] = deepcopy(value)
    return dict1


def update_nested_dict_node(dict_node: orm.Dict, update_dict: dict[str, Any], extend_list: bool = False) -> orm.Dict:
    """Utility to update a Dict node in a nested way"""
    pydict = dict_node.get_dict()
    update_nested_dict(pydict, update_dict, extend_list=extend_list)
    # Check if we have updated the node. If not, return the same node.
    if pydict == dict_node.get_dict():
        return dict_node
    return orm.Dict(dict=pydict)
