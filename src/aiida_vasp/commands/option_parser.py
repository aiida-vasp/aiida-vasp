"""
Module for reading options for cmd interface.

This module provides the :class:`OptionParser` class for parsing command-line options
in various formats including hierarchical dot notation, JSON, and YAML files.

**Key Features:**

* Hierarchical dot notation parsing (e.g., ``incar.elem=1,options.resources.num_machines=2``)
* Automatic type conversion (int, float, bool, None, str)
* JSON and YAML file loading support
* Backward compatibility functions

**Example Usage:**

.. code-block:: python

    from aiida_vasp.commands.option_parser import OptionParser

    # Parse hierarchical settings
    result = OptionParser.parse_hierarchical_dict("incar.elem=1,debug.enabled=true,solver.method=None")
    # Returns: {'incar': {'elem': 1}, 'debug': {'enabled': True}, 'solver': {'method': None}}

    # Process various option formats
    result = OptionParser.process_dict_option("config.json")  # Load from file
    result = OptionParser.process_dict_option('{"key": "value"}')  # Parse JSON
    result = OptionParser.process_dict_option("key=value,nested.key=123")  # Parse hierarchical
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import yaml


class OptionParser:
    """Parser for command-line options supporting various formats including hierarchical dot notation."""

    @staticmethod
    def nested_dict():
        """Create a nested defaultdict that automatically creates missing levels.

        :return: A nested defaultdict instance
        :rtype: defaultdict
        """
        return defaultdict(OptionParser.nested_dict)

    @classmethod
    def process_dict_option(cls, value: str | None) -> dict:
        """Process an option that can be a JSON string or a file path.

        :param value: String value that can be JSON, file path, or hierarchical key=value pairs
        :type value: str or None
        :return: Parsed dictionary from the input value
        :rtype: dict
        """
        if value is None:
            return {}
        if value.lower().endswith(('.json', '.yaml', '.yml')):
            return cls._load_dict_from_file(value)
        return cls._parse_text_as_dict(value)

    @classmethod
    def parse_hierarchical_dict(cls, settings_str: str) -> dict:
        """Parse hierarchical settings with dot notation into nested dictionaries.

        Uses defaultdict for more concise code and automatic creation of nested levels.

        :param settings_str: String containing comma-separated key=value pairs,
                            where keys can use dot notation for nesting
        :type settings_str: str
        :return: Nested dictionary structure based on the dot notation
        :rtype: dict

        **Examples:**

        .. code-block:: python

            >>> OptionParser.parse_hierarchical_dict(
            ...     "incar.elem=1,relax_settings.algo=rd,options.resources.num_machines=1"
            ... )
            {'incar': {'elem': 1}, 'relax_settings': {'algo': 'rd'}, 'options': {'resources': {'num_machines': 1}}}

            >>> OptionParser.parse_hierarchical_dict("debug.enabled=true,solver.tolerance=1e-6,mesh.nx=100")
            {'debug': {'enabled': True}, 'solver': {'tolerance': 1e-06}, 'mesh': {'nx': 100}}
        """
        if not settings_str:
            return {}

        result = cls.nested_dict()

        for pair in settings_str.split(','):
            if '=' not in pair:
                continue

            key_path, value = pair.split('=', 1)
            key_path = key_path.strip()
            value = value.strip()

            # Convert value to appropriate type
            converted_value = cls.convert_value(value)

            # Split the key path by dots and strip whitespace
            keys = [k.strip() for k in key_path.split('.')]

            # Navigate to the nested location using defaultdict
            current_dict = result
            for key in keys[:-1]:
                current_dict = current_dict[key]

            # Set the final value
            current_dict[keys[-1]] = converted_value

        # Convert defaultdict back to regular dict for clean output
        return cls._defaultdict_to_dict(result)

    @staticmethod
    def convert_value(value: str):
        """Convert string value to appropriate Python type.

        :param value: String value to convert
        :type value: str
        :return: Converted value (int, float, bool, None, or str)
        :rtype: int or float or bool or None or str

        **Supported conversions:**

        * Integers: ``"123"`` → ``123``
        * Floats: ``"1.5"`` → ``1.5``
        * Booleans: ``"true"``, ``"yes"``, ``"on"``, ``"1"`` → ``True``
        * Booleans: ``"false"``, ``"no"``, ``"off"``, ``"0"`` → ``False``
        * None: ``"None"``, ``"null"``, ``"nil"`` → ``None``
        * Strings: Everything else remains as string
        """
        if not value:
            return value

        # Try to convert to None
        if value.lower() in ('none', 'null', 'nil'):
            return None

        # Try to convert to int
        try:
            return int(value)
        except ValueError:
            pass

        # Try to convert to float
        try:
            return float(value)
        except ValueError:
            pass

        # Try to convert to boolean
        if value.lower() in ('true', 'yes', 'on', '1'):
            return True
        elif value.lower() in ('false', 'no', 'off', '0'):
            return False

        # Return as string
        return value

    @staticmethod
    def _defaultdict_to_dict(d):
        """Recursively convert defaultdict to regular dict.

        :param d: Input data structure to convert
        :return: Regular dictionary with all defaultdicts converted
        :rtype: dict or original type
        """
        if isinstance(d, defaultdict):
            return {k: OptionParser._defaultdict_to_dict(v) for k, v in d.items()}
        return d

    @classmethod
    def _parse_text_as_dict(cls, resources_str: str) -> dict:
        """Parse resources from various formats defined directly as a string.

        :param resources_str: String to parse (JSON or key=value format)
        :type resources_str: str
        :return: Parsed dictionary
        :rtype: dict
        """
        if not resources_str:
            return {}
        # Try JSON first
        try:
            return json.loads(resources_str)
        except json.JSONDecodeError:
            pass

        # Try key=value format
        return cls.parse_hierarchical_dict(resources_str)

    @staticmethod
    def _load_dict_from_file(overrides_path: Path | str) -> dict:
        """Load some settings from a file.

        :param overrides_path: Path to the file containing settings
        :type overrides_path: Path or str
        :return: Dictionary loaded from file
        :rtype: dict
        :raises: json.JSONDecodeError, yaml.YAMLError for malformed files
        """
        if not overrides_path:
            return {}

        overrides_path = Path(overrides_path)
        extension = overrides_path.suffix.lower()

        with open(overrides_path, 'r', encoding='utf8') as f:
            if extension in ['.json']:
                return json.load(f)
            elif extension in ['.yaml', '.yml']:
                return yaml.safe_load(f)
            else:
                # Try YAML first, then JSON
                content = f.read()
                try:
                    return yaml.safe_load(content)
                except yaml.YAMLError:
                    return json.loads(content)


# Convenience functions for backward compatibility
def process_dict_option(value: str | None) -> dict:
    """Backward compatibility wrapper for OptionParser.process_dict_option.

    :param value: String value that can be JSON, file path, or hierarchical key=value pairs
    :type value: str or None
    :return: Parsed dictionary from the input value
    :rtype: dict
    """
    return OptionParser.process_dict_option(value)


def parse_hierarchical_dict(settings_str: str) -> dict:
    """Backward compatibility wrapper for OptionParser.parse_hierarchical_dict.

    :param settings_str: String containing comma-separated key=value pairs
    :type settings_str: str
    :return: Nested dictionary structure based on the dot notation
    :rtype: dict
    """
    return OptionParser.parse_hierarchical_dict(settings_str)


def convert_value(value: str):
    """Backward compatibility wrapper for OptionParser.convert_value.

    :param value: String value to convert
    :type value: str
    :return: Converted value (int, float, bool, None, or str)
    :rtype: int or float or bool or None or str
    """
    return OptionParser.convert_value(value)
