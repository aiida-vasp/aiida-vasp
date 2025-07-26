"""
The .win parser interface.

=========================

Contains routines to parse Wannier90 input files. Will in the future utilize
the parser in the Wannier90 plugin, but no input parser exists yet.
"""

import re
from typing import Any, Callable


class BaseKeyValueParser:  # pylint: disable=useless-object-inheritance
    """
    Common codebase for all parser utilities.

    This class provides utility methods for parsing key-value and line-based data.
    """

    empty_line = re.compile(r'[\r\n]\s*[\r\n]')

    @classmethod
    def line(cls, fobj_or_str: str | Any, d_type: type = str) -> Any:
        """
        Grab a line from a file object or string and convert it to ``d_type`` (default: ``str``).

        :param fobj_or_str: File object or string to read the line from.
        :type fobj_or_str: file or str
        :param d_type: Type to convert each item in the line to. Defaults to ``str``.
        :type d_type: type
        :return: Single value or list of values, depending on the number of items in the line.
        """
        if isinstance(fobj_or_str, str):
            line = fobj_or_str
        else:
            line = fobj_or_str.readline()
        # previously this was map instead of list comprehension
        res = [d_type(item) for item in line.split()]
        if len(res) == 1:
            return res[0]
        return res

    @classmethod
    def splitlines(cls, fobj_or_str: str | Any, d_type: type = float) -> list[Any]:
        """
        Split a chunk of text into a list of lines and convert each line to ``d_type`` (default: ``float``).

        :param fobj_or_str: File object or string to split into lines.
        :type fobj_or_str: file or str
        :param d_type: Type to convert each item in the line to. Defaults to ``float``.
        :type d_type: type
        :return: List of values, one per line.
        """
        if isinstance(fobj_or_str, str):
            lines = fobj_or_str.split('\n')
        else:
            lines = fobj_or_str.readlines()
        return [cls.line(line, d_type) for line in lines]


class KeyValueParser(BaseKeyValueParser):
    """
    Key and value parser.

    This base class provides utility functions for parsing files that are
    (mostly) in a ``key = value`` format.

    .. note::
        This class does not integrate with the ``VaspParser`` class currently.

    Example usage::

        import re
        from aiida_vasp.parsers.file_parsers.parser import KeyValueParser

        class ParamParser(KeyValueParser):
            def __init__(self, file_path):
                self._file_path = py.path.local(file_path)
                super().__init__()
                self.result = {}
            def convert_or_not(self, value):
                for converter in self.get_converter_iter():
                    converted = self.try_convert(value, converter)
                    if converted and 'value' in converted:
                        return converted['value']
                return value
            def parse_file(self):
                assignments = re.findall(self.assignment, self._file_path.read())
                return {key: self.convert_or_not(value)}

    Parses files like::

        StrParam = value_1
        FloatParam = 1.0
        BoolParam = True
    """

    assignment = re.compile(r'(\w+)\s*[=: ]\s*([^;\n]*);?')
    bool_true = re.compile(r'^T$')
    bool_false = re.compile(r'^F$')
    comment = True

    @classmethod
    def get_lines(cls, filename: str) -> list[str]:
        """
        Read all lines from a file.

        :param filename: Path to the file.
        :type filename: str
        :return: List of lines from the file.
        :rtype: list[str]
        """
        with open(filename, 'r', encoding='utf8') as input_file:
            lines = input_file.read().splitlines()
        return lines

    @classmethod
    def retval(cls, *args: Any, **kwargs: Any) -> dict[str, Any]:
        """
        Normalize return values from value conversion functions.

        :return: Dictionary with the value and any additional keyword arguments.
        :rtype: dict
        """
        val = list(args)
        if len(val) == 1:
            val = val[0]
        ret = {'value': val}
        ret.update(kwargs)
        return ret

    @classmethod
    def flatten(cls, lst: list[list[Any]]) -> list[Any]:
        """
        Flatten a list of lists into a single list.

        :param lst: List of lists.
        :type lst: list
        :return: Flattened list.
        :rtype: list
        """
        return [i for j in lst for i in j]

    @classmethod
    def find_kv(cls, line: str) -> list[tuple[str, str]]:
        """
        Find key-value pairs in a line using the assignment regex.

        :param line: Line to search for key-value pairs.
        :type line: str
        :return: List of (key, value) tuples.
        :rtype: list[tuple]
        """
        return re.findall(cls.assignment, line)

    @classmethod
    def float(cls, string_: str) -> dict[str, Any]:
        """
        Parse a string into a float value followed by a comment.

        :param string_: String to parse.
        :type string_: str
        :return: Dictionary with value and comment.
        :rtype: dict
        """
        vals = string_.split()
        value = float(vals.pop(0))
        comment = ' '.join(vals)
        return cls.retval(value, comment=comment)

    @classmethod
    def float_unit(cls, string_: str) -> dict[str, Any]:
        """
        Parse string into a float number with attached unit.

        :param string_: String to parse.
        :type string_: str
        :return: Dictionary with value, unit, and comment.
        :rtype: dict
        """
        vals = string_.split()
        value = float(vals.pop(0))
        unit = vals.pop(0) if vals else ''
        comment = ' '.join(vals)
        return cls.retval(value, unit, comment=comment)

    @classmethod
    def int(cls, string_: str) -> dict[str, Any]:
        """
        Parse a string into an integer value followed by a comment.

        :param string_: String to parse.
        :type string_: str
        :return: Dictionary with value and comment.
        :rtype: dict
        """
        vals = string_.split()
        value = int(vals.pop(0))
        comment = ' '.join(vals)
        return cls.retval(value, comment=comment)

    @classmethod
    def int_unit(cls, string_: str) -> dict[str, Any]:
        """
        Convert a string into a python value, associated unit and optional comment.

        :param string_: String to parse.
        :type string_: str
        :return: Dictionary with value, unit, and comment.
        :rtype: dict
        """
        vals = string_.split()
        value = int(vals.pop(0))
        unit = vals.pop(0) if vals else ''
        comment = ' '.join(vals)
        return cls.retval(value, unit, comment=comment)

    @classmethod
    def string(cls, string_: str) -> dict[str, Any]:
        """
        Parse a string into value and comment, assuming only the first word is the value.

        :param string_: String to parse.
        :type string_: str
        :return: Dictionary with value and comment.
        :rtype: dict
        """
        vals = string_.split()
        value = vals.pop(0)
        comment = ' '.join(vals)
        return cls.retval(value, comment=comment)

    @classmethod
    def bool(cls, string_: str) -> dict[str, Any]:
        """
        Parse string into a boolean value.

        :param string_: String to parse.
        :type string_: str
        :return: Dictionary with value and comment.
        :rtype: dict
        :raises ValueError: If the string does not match a boolean pattern.
        """
        vals = string_.split()
        bool_str = vals.pop(0)
        if re.match(cls.bool_true, bool_str):
            value = True
        elif re.match(cls.bool_false, bool_str):
            value = False
        else:
            raise ValueError(
                f'bool string {string_} did not match any of {[cls.bool_true.pattern, cls.bool_false.pattern]}'
            )
        comment = ' '.join(vals)
        return cls.retval(value, comment=comment)

    @classmethod
    def kv_list(cls, filename: str) -> list[Any]:
        """
        Read a file and return a list of key-value pairs for each line.

        :param filename: Path to the file.
        :type filename: str
        :return: List of key-value pairs.
        :rtype: list
        """
        with open(filename, 'r', encoding='utf8') as input_fo:
            kv_list = filter(None, map(cls.find_kv, input_fo))
        return kv_list

    @classmethod
    def kv_dict(cls, kv_list: list[Any]) -> dict[str, Any]:
        """
        Convert a list of key-value pairs into a dictionary.

        :param kv_list: List of key-value pairs.
        :type kv_list: list
        :return: Dictionary of key-value pairs.
        :rtype: dict
        """
        kv_dict = dict(cls.flatten(kv_list))
        return kv_dict

    @classmethod
    def clean_value(cls, str_value: str) -> dict[str, Any]:
        """
        Get the converted python value from a string.

        :param str_value: String value to convert.
        :type str_value: str
        :return: Dictionary with the converted value.
        :rtype: dict
        """
        if str_value == '':
            return cls.retval(str_value)
        cleaned_value = None
        converters = cls.get_converter_iter()
        while not cleaned_value:
            cleaned_value = cls.try_convert(str_value, converters.next())
        return cleaned_value

    @classmethod
    def get_converter_iter(cls) -> Any:
        """
        Get an iterator over the value converter functions in order.

        :return: Iterator over converter functions.
        """
        converter_order = [cls.bool, cls.int, cls.float, cls.string]
        return (i for i in converter_order)

    @classmethod
    def try_convert(cls, input_value: str, converter: Callable[[str], dict[str, Any]]) -> dict[str, Any] | None:
        """
        Try to convert the input string into a python value given a conversion function.

        :param input_value: Value to convert.
        :type input_value: str
        :param converter: Converter function to use.
        :type converter: callable
        :return: Dictionary with the converted value, or None if conversion failed.
        :rtype: dict or None
        """
        if not isinstance(input_value, str):
            return {'value': input_value}
        try:
            cleaned_value = converter(input_value)
        except ValueError:
            cleaned_value = {}

        if cleaned_value.get('value', None) is None:
            return None
        return cleaned_value


class WinParser(KeyValueParser):
    """
    Parses wannier90 input files.

    This parser extracts keywords, blocks, and comments from a Wannier90 ``.win`` input file.
    """

    block = re.compile(r'begin (?P<name>\w*)\s*\n\s*(?P<content>[\w\W]*)\s*\n\s*end \1')
    comment = re.compile(r'(!.*)\n?')

    def __init__(self, path: str) -> None:
        """
        Initialize the parser and parse the Wannier90 input file.

        :param path: Path to the Wannier90 .win file.
        :type path: str
        """
        super().__init__()
        self.result = {}
        with open(path, 'r', encoding='utf8') as winf:
            self.keywords, self.blocks, self.comments = WinParser.parse_win(winf)
        self.result.update(self.keywords)
        self.result.update(self.blocks)

    @classmethod
    def parse_win(cls, fobj_or_str: str | Any) -> tuple[dict[str, Any], dict[str, list[str]], list[str]]:
        """
        Parse a Wannier90 input file or string.

        :param fobj_or_str: File object or string containing the Wannier90 input.
        :type fobj_or_str: file or str
        :return: Tuple of (keywords dict, blocks dict, comments list).
        :rtype: tuple
        """
        if isinstance(fobj_or_str, str):
            content = fobj_or_str
        else:
            content = fobj_or_str.read()
        comments = re.findall(cls.comment, content)
        content = re.sub(cls.comment, '', content)
        blocks = re.findall(cls.block, content)
        content = re.sub(cls.block, '', content)
        kvd = dict(re.findall(cls.assignment, content))
        bld = {}
        for keyword, value in blocks:
            # do not split individual lines
            bld[keyword] = [line.strip() for line in value.split('\n')]
        return kvd, bld, comments
