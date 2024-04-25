"""
The .win parser interface.

--------------------------
Contains routines to parse Wannier90 input files. Will in the future be utilizing
the parser in the Wannier90 plugin, but no input parser exists.
"""

import re


class BaseKeyValueParser:  # pylint: disable=useless-object-inheritance
    """Common codebase for all parser utilities."""

    empty_line = re.compile(r'[\r\n]\s*[\r\n]')

    @classmethod
    def line(cls, fobj_or_str, d_type=str):
        """Grab a line from a file object or string and convert it to d_type (default: str)."""
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
    def splitlines(cls, fobj_or_str, d_type=float):
        """Split a chunk of text into a list of lines and convert each line to d_type (default: float)."""
        if isinstance(fobj_or_str, str):
            lines = fobj_or_str.split('\n')
        else:
            lines = fobj_or_str.readlines()
        return [cls.line(line, d_type) for line in lines]


class KeyValueParser(BaseKeyValueParser):
    """
    Key and value parser.
    ---------------------
    This baseclass has some utility functions for parsing files that are
    (mostly) in a `key` = `value` format.
    This class does not integrate with the VaspParser class currently.
    A simple example, which tries to convert values to python objects on a best effort basis. Usage::
        import re
        from aiida_vasp.parsers.file_parsers.parser import KeyValueParser
        ParamParser(KeyValueParser):
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
    def get_lines(cls, filename):
        with open(filename, 'r', encoding='utf8') as input_file:
            lines = input_file.read().splitlines()
        return lines

    @classmethod
    def retval(cls, *args, **kwargs):
        """Normalize return values from value conversion functions."""
        val = list(args)
        if len(val) == 1:
            val = val[0]
        ret = {'value': val}
        ret.update(kwargs)
        return ret

    @classmethod
    def flatten(cls, lst):
        return [i for j in lst for i in j]

    @classmethod
    def find_kv(cls, line):
        return re.findall(cls.assignment, line)

    @classmethod
    def float(cls, string_):
        """Parse a string into an float value followed by a comment."""
        vals = string_.split()
        value = float(vals.pop(0))
        comment = ' '.join(vals)
        return cls.retval(value, comment=comment)

    @classmethod
    def float_unit(cls, string_):
        """Parse string into a float number with attached unit."""
        vals = string_.split()
        value = float(vals.pop(0))
        unit = vals.pop(0) if vals else ''
        comment = ' '.join(vals)
        return cls.retval(value, unit, comment=comment)

    @classmethod
    def int(cls, string_):
        """Parse a string into an integer value followed by a comment."""
        vals = string_.split()
        value = int(vals.pop(0))
        comment = ' '.join(vals)
        return cls.retval(value, comment=comment)

    @classmethod
    def int_unit(cls, string_):
        """Convert a string into a python value, associated unit and optional comment."""
        vals = string_.split()
        value = int(vals.pop(0))
        unit = vals.pop(0) if vals else ''
        comment = ' '.join(vals)
        return cls.retval(value, unit, comment=comment)

    @classmethod
    def string(cls, string_):
        """Parse a string into value and comment, assuming only the first word is the value."""
        vals = string_.split()
        value = vals.pop(0)
        comment = ' '.join(vals)
        return cls.retval(value, comment=comment)

    @classmethod
    def bool(cls, string_):
        """Parse string into a boolean value."""
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
    def kv_list(cls, filename):
        with open(filename, 'r', encoding='utf8') as input_fo:
            kv_list = filter(None, map(cls.find_kv, input_fo))
        return kv_list

    @classmethod
    def kv_dict(cls, kv_list):
        kv_dict = dict(cls.flatten(kv_list))
        return kv_dict

    @classmethod
    def clean_value(cls, str_value):
        """Get the converted python value from a string."""
        if str_value == '':
            return cls.retval(str_value)
        cleaned_value = None
        converters = cls.get_converter_iter()
        while not cleaned_value:
            cleaned_value = cls.try_convert(str_value, converters.next())
        return cleaned_value

    @classmethod
    def get_converter_iter(cls):
        converter_order = [cls.bool, cls.int, cls.float, cls.string]
        return (i for i in converter_order)

    @classmethod
    def try_convert(cls, input_value, converter):
        """Try to convert the input string into a python value given a conversion function."""
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
    """Parses wannier90 input."""

    block = re.compile(r'begin (?P<name>\w*)\s*\n\s*(?P<content>[\w\W]*)\s*\n\s*end \1')
    comment = re.compile(r'(!.*)\n?')

    def __init__(self, path):  # pylint: disable=missing-function-docstring
        super().__init__()
        self.result = {}
        with open(path, 'r', encoding='utf8') as winf:
            self.keywords, self.blocks, self.comments = WinParser.parse_win(winf)
        self.result.update(self.keywords)
        self.result.update(self.blocks)

    @classmethod
    def parse_win(cls, fobj_or_str):
        """Parse a wannier90 input."""
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
