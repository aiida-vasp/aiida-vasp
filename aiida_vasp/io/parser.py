"""
This module contains the base class for other vasp parsers.
"""
import re

from six import string_types


class BaseParser(object):
    """Common codebase for all parser utilities. It provides the following interface to be used
    by the VaspParser:

        - properties: a list holding all the properties this parser can extract from it's file
        - get_quantities(properties, output): Method to be called by the VaspParser
          which will call the specific parsing methods for each requested property and add
          their results to the correct ouput node.

          :output contains data parsed by other file parsers and optionally a 'settings' card
                  which determines the behaviour of each file parsers _parse_file method.

    """
    empty_line = re.compile(r'[\r\n]\s*[\r\n]')

    def __init__(self):
        super(BaseParser, self).__init__()
        self._parsable_items = {}
        self._parsed_data = None
        self._filename = None
        self._filepath = None
        self._out_folder = None

    @classmethod
    def line(cls, fobj_or_str, d_type=str):
        """Grab a line from a file object or string and convert it to d_type (default: str)"""
        if isinstance(fobj_or_str, str):
            line = fobj_or_str
        else:
            line = fobj_or_str.readline()
        res = map(d_type, line.split())
        if len(res) == 1:
            return res[0]
        return res

    @classmethod
    def splitlines(cls, fobj_or_str, d_type=float):
        """split a chunk of text into a list of lines and convert each line to d_type (default: float)"""
        if isinstance(fobj_or_str, str):
            lines = fobj_or_str.split('\n')
        else:
            lines = fobj_or_str.readlines()
        return [cls.line(l, d_type) for l in lines]

    def get_quantity(self, quantity, output):
        """
        Public method to get the required quantity from the _parsed_data dictionary if that exists.
        Otherwise parse the file. This method will be registered to the VaspParsers get_quantities
        delegate during __init__.
        """
        if quantity in self._parsable_items:
            # gather everything required for parsing this component
            inputs = {}
            inputs['settings'] = output.get('settings')
            for inp in self._parsable_items[quantity]['inputs']:
                inputs[inp] = output.get(inp)

            if not self._parsed_data:
                # The file has not been parsed yet.
                self._parsed_data = self._parse_file(inputs)

            output[quantity] = self._parsed_data.get(quantity)

    def _parse_file(self, inputs):
        """
        Abstract base method to parse this file parsers file. Has to be overwritten by the child class.
        """

        raise NotImplementedError('{0} does not implement a _parse_file() method.'.format(self.__class__.__name__))


# pylint: disable=abstract-method
class KeyValueParser(BaseParser):
    """
    contains regex and functions to find grammar elements
    in vasp input and output files
    """
    assignment = re.compile(r'(\w+)\s*[=: ]\s*([^;\n]*);?')
    bool_true = re.compile(r'^T$')
    bool_false = re.compile(r'^F$')
    comments = True

    @classmethod
    def get_lines(cls, filename):
        with open(filename) as input_file:
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
        """Parse string into a float number with attached unit"""
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
        """Parse string into a boolean value"""
        vals = string_.split()
        bool_str = vals.pop(0)
        if re.match(cls.bool_true, bool_str):
            value = True
        elif re.match(cls.bool_false, bool_str):
            value = False
        else:
            raise ValueError('bool string {} did not match any of {}'.format(string_, [cls.bool_true.pattern, cls.bool_false.pattern]))
        comment = ' '.join(vals)
        return cls.retval(value, comment=comment)

    @classmethod
    def kv_list(cls, filename):
        with open(filename) as input_fo:
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
        if not isinstance(input_value, string_types):
            return {'value': input_value}
        try:
            cleaned_value = converter(input_value)
        except ValueError:
            cleaned_value = {}

        if cleaned_value.get('value', None) is None:
            return None
        return cleaned_value
