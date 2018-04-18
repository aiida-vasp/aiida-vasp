r"""
Base classes for vasp file parsers.

BaseFileParser
--------------
This baseclass makes no assumptions about the format of the file to be parsed.

Usage::

    import re

    import py

    from aiida_vasp.utils.aiida_utils import get_data_class
    from aiida_vasp.io.parser import BaseFileParser

    ExampleFileParser(BaseFileParser):

        PARSABLE_ITEMS = {
            'item1': {
                'inputs': ['required_quantity'],  # This quantity will be parsed first and made available in time if possible
                'parsers': ['ExampleFile'], # During setup the VaspParser will check, whether ExampleFile has been retrieved
                ## and initialise the corresponding parser, if this quantity is requested by setting any of the
                ## 'parser_settings['add_OutputNode'] = True'
                'nodeName': ['examples'],  # The quantity will be added to the 'output_examples' output node
                'prerequisites: ['required_quantity'],  # This prohibits the parser from trying to parse item1 without ``required_quantity``
            }
            'item2': {
                'inputs': [],
                'parsers': ['ExampleFile'],
                'nodeName': ['examples'],
                'prerequisites': [],
            }
        }

        def __init__(self, *args, **kwargs):
            super(ExampleFileParser, self).__init__(*args, **kwargs)
            self._parsable_items = self.PARSABLE_ITEMS
            self._parsed_data = {}

        def _parse_file(self, inputs):
            example_file = py.path.local(self._file_path)  # self._file_path is set by the superclass
            data = example_file.read()
            item1 = int(re.findall(r'item1 is: (\d+)', data)[0]) * inputs['required_quantity']  # extract item 1
            item2 = [int(i) for i in re.findall(r'item2: (\d+)', data)]  # extract list of item2
            output_node = get_data_class('parameter')(dict={  # construct ParameterData node
                'item1': item1,
                'item2': item2
            }
            return {'examples': output_node}  # each of the ``nodeName``s from above must be a key in the returned dict

    example_parser = ExampleFileParser('example_file')
    item1 = example_parser.get_quantity('item1', inputs={'required_quantity': 1})
    item2 = example_parser.get_quantities('item2, inputs=None)

Parses Files like::

    item1 is: 213
    item2: 1
    item2: 2
    item2: 4

KeyValueParser
--------------
This baseclass has some utility functions for parsing files that are (mostly) in a `key` = `value` format.

This class does not integrate with the VaspParser class currently.

A simple example, which tries to convert values to python objects on a best effort basis:
Usage::

    import re

    import py

    from aiida_vasp.io.parser import KeyValueParser

    ParamParser(KeyValueParser):

        def __init__(self, file_path):
            self._file_path = py.path.local(file_path)
            super(WinParser, self).__init__()
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
import re

from six import string_types
from aiida_vasp.utils.delegates import delegate_method_kwargs


class BaseParser(object):
    """Common codebase for all parser utilities"""
    empty_line = re.compile(r'[\r\n]\s*[\r\n]')

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
        """Split a chunk of text into a list of lines and convert each line to d_type (default: float)."""
        if isinstance(fobj_or_str, str):
            lines = fobj_or_str.split('\n')
        else:
            lines = fobj_or_str.readlines()
        return [cls.line(l, d_type) for l in lines]


class BaseFileParser(BaseParser):
    """
    Abstract base class for the individual file parsers. It provides the following interface to be used by the VaspParser:

        - _parsable_items: a dictionary holding all items this parser can extract from it's file as well
          as the required information on how to extract those.
        - _parsed_data: a dictionary containing all the parsed data from this file.
        - get_quantity(): Method to be called by the VaspParser
          which will either fill the _parsed_data in case that it is empty by calling _parse_file
          or return the requested data from the _parsed_data. If another quantity is required as
          prerequisite it will be requested from the VaspParser.

          This method will be subscribed to the VaspParsers get_quantity delegate during initialisation.
          When the VaspParser calls his delegate this method will be called and return the requested
          quantity.
        _ _parse_file: an abstract method to be implemented by the actual file parser, which will
          parse the file and fill the _parsed_data dictionary.

        :param calc_parser_cls: Python class, optional, class of the calling CalculationParser instance

    The second way to use the BaseFileParser is for writing VASP files based on given Aiida data objects.
    The BaseFileParser will store the data object in self._data_obj and provide a public 'write' method
    to right the corresponding VASP file.

    For now the BaseFileParser has four different representations of it's data:
    - _filepath: data stored in a file
    . _parsed_data: dictionary with everything parsed from the file.
    - _data_obj: The Aiida data obj, which should be written to a file.
    - _parsed_obj: The parsevasp object representing this file.

    It seams like _parsed_data and _parsed_obj mgiht be equivalent. The same is true for _file_path and
    data_obj since both of them represent the data in its unparsed state. A future rework might unify this.
    also 'parse_file' could then switch it's parsing direction based on the context.
    """

    def __init__(self, file_path=None, calc_parser_cls=None):
        super(BaseFileParser, self).__init__()
        self._vasp_parser = calc_parser_cls
        if calc_parser_cls is not None:
            calc_parser_cls.get_quantity.add_listener(self.get_quantity)

        self._parsable_items = {}
        self._parsed_data = {}
        self._filepath = None
        self._data_obj = None
        self._parsed_obj = None

    @delegate_method_kwargs(prefix='_init_with_')
    def init_with_kwargs(self, **kwargs):
        """Delegate initialization to _init_with - methods."""

    def get_quantity(self, quantity, settings, inputs=None):
        """
        Public method to get the required quantity from the _parsed_data dictionary if that exists.

        Otherwise parse the file. This method will be registered to the VaspParsers get_quantities
        delegate during __init__.
        """

        if quantity not in self._parsable_items:
            return None

        if self._parsed_data is None or self._parsed_data.get(quantity) is None:
            # The file has not been parsed yet, or the quantity has not
            # been parsed yet, due to lack of required inputs..

            # gather everything required for parsing this component.
            if inputs is None:
                inputs = {}
            inputs['settings'] = settings

            if self._vasp_parser is not None:
                # gather everything required for parsing this quantity from the VaspParser.
                for inp in self._parsable_items[quantity]['inputs']:
                    inputs[inp] = self._vasp_parser.get_inputs(inp)
                    if inputs[inp] is None and inp in self._parsable_items[quantity]['prerequisites']:
                        # The VaspParser was unable to provide the required input.
                        return {quantity: None}

            self._parsed_data = self._parse_file(inputs)

        return {quantity: self._parsed_data.get(quantity)}

    def write(self, filepath):
        if self._data_obj is not None:
            self._data_obj.write(filepath)

    def _parse_file(self, inputs):
        """Abstract base method to parse this file parsers file. Has to be overwritten by the child class."""

        raise NotImplementedError('{0} does not implement a _parse_file() method.'.format(self.__class__.__name__))


class KeyValueParser(BaseParser):
    """Contains regex and functions to find grammar elements in vasp input and output files."""
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
