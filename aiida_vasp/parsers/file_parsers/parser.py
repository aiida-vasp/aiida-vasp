"""
Base classes for the VASP file parsers.

---------------------------------------
Contains the base classes for the VASP file parsers.
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


class BaseFileParser(BaseParser):
    """
    Abstract base class for the individual file parsers.

    It provides the following interface to be used by the VaspParser:

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

        :keyword file_path: Initialise with a path to a file. The file will be parsed by the FileParser
        :keyword data: Initialise with an aiida data object. This may be SingleFileData, KpointsData or StructureData.

    Additional keyword arguments might be defined by the inheriting classes.

    The second way to use the BaseFileParser is for writing VASP files based on given Aiida data objects.
    The BaseFileParser will provide the data object by the _parsed_object property and offer a public 'write' method
    to write the corresponding VASP file. Depending on whether the file under consideration is an actual input
    file, this may simply mean copying a file.
    """

    PARSABLE_ITEMS = {}

    def __init__(self, calc_parser_cls=None, **kwargs):  # pylint: disable=unused-argument
        super(BaseFileParser, self).__init__()
        self._logger = aiidalogger.getChild(self.__class__.__name__)
        self._vasp_parser = calc_parser_cls
        self.settings = None
        if calc_parser_cls is not None:
            calc_parser_cls.get_quantity.append(self.get_quantity)
            self.settings = calc_parser_cls.settings

        self.parsable_items = {}
        self._parsed_data = {}
        self._data_obj = None

    @delegate_method_kwargs(prefix='_init_with_')
    def init_with_kwargs(self, **kwargs):
        """Delegate initialization to _init_with - methods."""

    def _init_with_file_path(self, path):
        """Init with a file path."""
        self._data_obj = SingleFile(path=path)
        self.parsable_items = self.__class__.PARSABLE_ITEMS
        self._parsed_data = {}

    def _init_with_data(self, data):
        """
        Init with aiida-data.

        This has to be overriden by every FileParser, which deals with an
        Aiida data class other than SingleFileData.
        """

        self._data_obj = SingleFile(data=data)
        self.parsable_items = self.__class__.PARSABLE_ITEMS
        self._parsed_data = {}

    def _init_with_settings(self, settings):
        """Init with settings."""
        self.settings = settings

    def get_quantity(self, quantity, inputs=None):
        """
        Public method to get the required quantity from the _parsed_data dictionary if that exists.

        Otherwise parse the file. This method will be registered to the VaspParsers get_quantities
        delegate during __init__.
        """

        if quantity not in self.parsable_items:
            return None

        if self._parsed_data is None or self._parsed_data.get(quantity) is None:
            # The file has not been parsed yet, or the quantity has not
            # been parsed yet, due to lack of required inputs..

            # gather everything required for parsing this component.
            if inputs is None:
                inputs = {}

            if self._vasp_parser is not None:
                # gather everything required for parsing this quantity from the VaspParser.
                for inp in self.parsable_items[quantity]['inputs']:
                    inputs.update(self._vasp_parser.get_inputs(inp))
                    if inputs[inp] is None and inp in self.parsable_items[quantity]['prerequisites']:
                        # The VaspParser was unable to provide the required input.
                        return {quantity: None}
            self._parsed_data = self._parse_file(inputs)

        return {quantity: self._parsed_data.get(quantity)}

    def write(self, file_path):
        """
        Writes a VASP style file from the parsed Object.

        For non input files this means simply copying the file.
        """
        if self._parsed_object is not None:
            self._parsed_object.write(file_path)

    @property
    def _parsed_object(self):
        """
        Property to return the FileParsers _data_obj.

        The data_obj is either an instance of one of the parsevasp parser classes,
        which provide a write function, an instance of an aiida data node or an
        instance of SingleFile in case that it is just a file and does not have
        it's own 'write' method.

        In particular FileParsers storing aiida data nodes. will have to override this.
        """
        return self._data_obj

    @property
    def data_obj(self):
        return self._data_obj

    def _parse_file(self, inputs):
        """Abstract base method to parse this file parsers file. Has to be overwritten by the child class."""

        raise NotImplementedError('{classname} does not implement a _parse_file() ' 'method.'.format(classname=self.__class__.__name__))


class SingleFile(object):  # pylint: disable=useless-object-inheritance
    """
    Datastructure for a singleFile file providing a write method.

    This should get replaced, as soon as parsevasp has a dedicated class.
    """

    def __init__(self, **kwargs):
        super(SingleFile, self).__init__()
        self._path = None
        self._data = None
        self.init_with_kwargs(**kwargs)

    @delegate_method_kwargs(prefix='_init_with_')
    def init_with_kwargs(self, **kwqargs):
        """Delegate initialization to _init_with - methods."""

    def _init_with_path(self, path):
        self._path = path

    def _init_with_data(self, data):
        """Initialise with SingleFileData."""
        self._data = data

    @property
    def path(self):
        return self._path

    def write(self, dst):
        """Copy file to destination."""
        if self._path is not None:
            import shutil
            shutil.copyfile(self._path, dst)
            return

        if self._data is not None:
            with self._data.open() as input_obj, open(dst) as output_obj:
                lines = input_obj.readlines()
                output_obj.writelines(lines)


class KeyValueParser(BaseParser):
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
        if not isinstance(input_value, str):
            return {'value': input_value}
        try:
            cleaned_value = converter(input_value)
        except ValueError:
            cleaned_value = {}

        if cleaned_value.get('value', None) is None:
            return None
        return cleaned_value
