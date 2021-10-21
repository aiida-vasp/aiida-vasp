"""
Base classes for the VASP object parsers.

---------------------------------------
Contains the base classes for the VASP object parsers.
"""
# pylint: disable=import-outside-toplevel
import re
from aiida.common import AIIDA_LOGGER as aiidalogger
from aiida_vasp.utils.delegates import delegate_method_kwargs


class BaseParser(object):  # pylint: disable=useless-object-inheritance
    """Common codebase for all parser utilities."""
    empty_line = re.compile(r'[\r\n]\s*[\r\n]')

    @classmethod
    def line(cls, fobj_or_str, d_type=str):
        """Grab a line from an object or string and convert it to d_type (default: str)."""
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


class BaseFileParser(BaseParser):
    """
    Abstract base class for the individual object parsers.

    It provides the following interface to be used by the VaspParser:

        - _parsable_items: a dictionary holding all items this parser can extract from it's object as well
          as the required information on how to extract those.
        - _parsed_data: a dictionary containing all the parsed data from this object.
        - get_quantity(): Method to be called by the VaspParser
          which will either fill the _parsed_data in case that it is empty by calling _parse_object
          or return the requested data from the _parsed_data. If another quantity is required as
          prerequisite it will be requested from the VaspParser.

          This method will be subscribed to the VaspParsers get_quantity delegate during initialisation.
          When the VaspParser calls his delegate this method will be called and return the requested
          quantity.
        _ _parse_object: an abstract method to be implemented by the actual object parser, which will
          parse the object and fill the _parsed_data dictionary.

        :param calc_parser_cls: Python class, optional, class of the calling CalculationParser instance

        :keyword path: Initialise with a path to a object. The object will be parsed by the object parser
        :keyword data: Initialise with an aiida data object. This may be SingleFileData, KpointsData or StructureData.

    Additional keyword arguments might be defined by the inheriting classes.

    The second way to use the BaseFileParser is for writing VASP objects based on given Aiida data objects.
    The BaseFileParser will provide the data object by the _parsed_object property and offer a public 'write' method
    to write the corresponding VASP object. Depending on whether the object under consideration is an actual input
    object, this may simply mean copying it.
    """

    PARSABLE_ITEMS = {}

    def __init__(self, **kwargs):  # pylint: disable=unused-argument
        super(BaseFileParser, self).__init__()
        self._logger = aiidalogger.getChild(self.__class__.__name__)
        self._exit_code = None
        self._parsable_items = self.PARSABLE_ITEMS
        self._parsed_data = {}
        if 'handler' in kwargs:
            # For instance a file handler
            self._data_obj = SingleFile(handler=kwargs['handler'])
        elif 'data' in kwargs:
            # An AiiDA data structure
            self._data_obj = SingleFile(data=kwargs['data'])
        else:
            self._data_obj = None

    @property
    def parsable_items(self):
        return self._parsable_items

    @property
    def exit_code(self):
        return self._exit_code

    def get_quantity(self, quantity_key):
        """
        Get the required quantity

        Either from the _parsed_data dictionary if that exists, otherwise parse the object.
        """

        if quantity_key not in self._parsable_items:
            return None

        if self._parsed_data.get(quantity_key) is None:
            self._parsed_data = self._parse_object({})

        return self._parsed_data.get(quantity_key)

    def get_quantity_from_inputs(self, quantity_name, inputs, vasp_parser):
        """Method to handle inputs (to be removed)"""

        if quantity_name not in self._parsable_items:
            return None

        if self._parsed_data.get(quantity_name) is None:
            if inputs is None:
                inputs = {}
            if vasp_parser is not None:
                # gather everything required for parsing this quantity from the VaspParser.
                for inp in self._parsable_items[quantity_name]['inputs']:
                    inputs.update(vasp_parser.get_inputs(inp))
                    if inputs[inp] is None and inp in self._parsable_items[quantity_name]['prerequisites']:
                        # The VaspParser was unable to provide the required input.
                        return None
            self._parsed_data = self._parse_object(inputs)

        return self._parsed_data.get(quantity_name)

    def write(self, path):
        """
        Writes a VASP object from the parsed object.

        For non input objects this means simply copying the object.
        """
        if self._parsed_object is not None:
            self._parsed_object.write(path)

    @property
    def _parsed_object(self):
        """
        Property to return the object parsers _data_obj.

        The data_obj is either an instance of one of the parsevasp parser classes,
        which provide a write function, an instance of an AiiDA data node or an
        instance of SingleFile in case it is just an object and does not have
        it's own 'write' method.

        In particular the specific object parsers storing AiiDA data nodes have to override this.
        """
        return self._data_obj

    @property
    def data_obj(self):
        return self._data_obj

    def _parse_object(self, inputs):
        """Abstract base method to parse this object. Has to be overwritten by the child class."""

        raise NotImplementedError('{classname} does not implement a _parse_object() ' 'method.'.format(classname=self.__class__.__name__))


class SingleFile(object):  # pylint: disable=useless-object-inheritance
    """
    Datastructure for a SingleFile object providing a write method.

    This should get replaced, as soon as parsevasp has a dedicated class.
    """

    def __init__(self, **kwargs):
        super(SingleFile, self).__init__()
        self._handler = None
        self._data = None
        self.init_with_kwargs(**kwargs)

    @delegate_method_kwargs(prefix='_init_with_')
    def init_with_kwargs(self, **kwqargs):
        """Delegate initialization to _init_with - methods."""

    def _init_with_handler(self, handler):
        self._handler = handler

    def _init_with_data(self, data):
        """Initialise with SingleFileData."""
        self._data = data

    @property
    def handler(self):
        return self._handler

    def write(self, dst):
        """Copy to destination."""
        if self._handler is not None:
            with open(dst, 'w') as output_obj:
                lines = self._handler.readlines()
                output_obj.writelines(lines)
            return

        if self._data is not None:
            with self._data.open() as input_obj, open(dst) as output_obj:
                lines = input_obj.readlines()
                output_obj.writelines(lines)


class KeyValueParser(BaseParser):
    """
    Key and value parser.

    ---------------------
    This baseclass has some utility functions for parsing objects that are
    (mostly) in a `key` = `value` format.

    This class does not integrate with the VaspParser class currently.

    A simple example, which tries to convert values to python objects on a best effort basis. Usage::

        import re

        from aiida_vasp.parsers.content_parsers.parser import KeyValueParser

        ParamParser(KeyValueParser):

            def __init__(self, path):
                self._path = py.path.local(path)
                super(WinParser, self).__init__()
                self.result = {}

            def convert_or_not(self, value):
                for converter in self.get_converter_iter():
                    converted = self.try_convert(value, converter)
                    if converted and 'value' in converted:
                        return converted['value']
                return value

            def parse_object(self):
                assignments = re.findall(self.assignment, self._path.read())
                return {key: self.convert_or_not(value)}

    Parses objects like::

        StrParam = value_1
        FloatParam = 1.0
        BoolParam = True
    """
    assignment = re.compile(r'(\w+)\s*[=: ]\s*([^;\n]*);?')
    bool_true = re.compile(r'^T$')
    bool_false = re.compile(r'^F$')
    comments = True

    @classmethod
    def get_lines(cls, name):
        with open(name) as input_object:
            lines = input_object.read().splitlines()
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
            raise ValueError('bool string {} did not match any of {}'.format(string_, [cls.bool_true.pattern, cls.bool_false.pattern]))
        comment = ' '.join(vals)
        return cls.retval(value, comment=comment)

    @classmethod
    def kv_list(cls, name):
        with open(name) as input_fo:
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
