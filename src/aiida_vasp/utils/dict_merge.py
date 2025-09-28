"""
Functions to merge dictionaries
"""

import collections
from copy import deepcopy


def recursive_merge_orig(left: dict, right: dict) -> dict:
    """Recursively merge two dictionaries into a single dictionary.

    If any key is present in both ``left`` and ``right`` dictionaries, the value from the ``right`` dictionary is
    assigned to the key.

    :param left: first dictionary
    :param right: second dictionary
    :return: the recursively merged dictionary
    """

    # Note that a deepcopy is not necessary, since this function is called recusively.
    right = right.copy()

    for key, value in left.items():
        if key in right:
            if isinstance(value, collections.abc.Mapping) and isinstance(right[key], collections.abc.Mapping):
                right[key] = recursive_merge_orig(value, right[key])

    merged = left.copy()
    merged.update(right)

    return merged


def recursive_merge(left: dict, right: dict) -> dict:
    """
    Recursively merge two dictionaries into a single dictionary, supporting special operations for keys.

    If a key is present in both ``left`` and ``right`` dictionaries, the value from ``right`` is used.
    If both values are dictionaries, they are merged recursively.

    Special operations (when the value in ``right`` is a dictionary with one of these keys):

    - ``$!del``: Delete this key from the result.
    - ``$!append``: Append a value to a list at this key.
    - ``$!extend``: Extend a list at this key with another list.
    - ``$!replace``: Replace the value at this key entirely.

    :param dict left: First dictionary.
    :param dict right: Second dictionary.
    :return: The recursively merged dictionary.
    :rtype: dict

    **Examples**

    Basic merge::

        >>> recursive_merge({'a': 1, 'b': {'c': 2}}, {'b': {'c': 3}, 'd': 4})
        {'a': 1, 'b': {'c': 3}, 'd': 4}

    Nested merge and list extend::

        >>> recursive_merge({'a': {'x': 1, 'y': {'z': 2}}, 'b': [1, 2]}, {'a': {'y': {'z': 3}}, 'b':
        {'$!extend': [3, 4]}})
        {'a': {'x': 1, 'y': {'z': 3}}, 'b': [1, 2, 3, 4]}

    Special: append and delete::

        >>> recursive_merge({'a': [1, 2], 'b': {'c': 5}}, {'a': {'$!append': 3}, 'b': {'$!del': True}})
        {'a': [1, 2, 3]}

        >>> recursive_merge({'a': [1, 2], 'b': {'c': 5}}, {'a': {'$!append': 3}, 'b': '$!del'}})
        {'a': [1, 2, 3]}

    Replace::

        >>> recursive_merge({'a': {'x': 1}}, {'a': {'$!replace': {'y': 2}}})
        {'a': {'y': 2}}
    """
    # Here a deepcopy is necessary as in-place modification is used
    left = deepcopy(left)
    for key, value_right in right.items():
        if key in left:
            # Apply special operations
            if '$!del' == value_right:
                del left[key]
                continue
            if isinstance(value_right, dict):
                if '$!del' in value_right:
                    del left[key]
                    continue
                if '$!append' in value_right:
                    left[key] = left[key].copy()
                    left[key].append(value_right['$!append'])
                    continue
                if '$!extend' in value_right:
                    left[key] = left[key].copy()
                    left[key].extend(value_right['$!extend'])
                    continue
                if '$!replace' in value_right:
                    left[key] = value_right['$!replace']
                    continue
            # Nested update
            if isinstance(value_right, collections.abc.Mapping) and isinstance(left[key], collections.abc.Mapping):
                left[key] = recursive_merge(left[key], value_right)
            else:
                left[key] = value_right
        else:
            left[key] = value_right
    return left
