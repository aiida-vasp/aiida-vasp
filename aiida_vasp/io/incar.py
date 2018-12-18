"""Provides INCAR file interface and utilities."""
import re
import operator
import functools
from collections import OrderedDict
import numpy as np
import pyparsing as pp

import six
from pymatgen.io.vasp import Incar as IncarPymatgen
from parsevasp.incar import Incar as IncarParsevasp
from aiida_vasp.io.parser import BaseFileParser

from aiida_vasp.utils.aiida_utils import get_data_node, get_data_class


class IncarItem(object):
    """
    Represent an item in an incar file.

    Example::

        ENCUT = 350 eV this is an example
        ^       ^   ^-------------------^
        name    value       comment

    :param name: mandatory
    :param value: mandatory
    :param comment: optional
    """
    _STR_TPL = '{name} = {value}{comment_sep}{comment}'

    def __init__(self, *args, **kwargs):
        if args and isinstance(args[0], self.__class__):
            self.set_copy(args[0])
        elif 'incar_item' in kwargs:
            self.set_copy(kwargs['incar_item'])
        else:
            self.value = None
            self._name = None
            self._comment = None
            self.set_values(*args, **kwargs)

    @classmethod
    def from_string(cls, input_string):
        result = IncarParamParser.get_parser().parseString(input_string)
        return IncarItem(result.name, result.value, result.comment)

    def set_copy(self, other):
        self.set_values(other.name, other.value, other.comment)

    def set_values(self, name, value, comment=None):
        self.name = name
        self.value = value
        self.comment = comment

    @property
    def comment(self):
        return self._comment if self._comment else ''

    @comment.setter
    def comment(self, comment_input):
        if comment_input is not None:
            self._comment = str(comment_input).lstrip('#').strip()

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name_input):
        self._name = name_input.upper()

    @property
    def comment_separator(self):
        if self.comment:
            return ' # '
        return ''

    def as_dict(self):
        return {'name': self.name, 'value': self.value, 'comment': self.comment}

    def __str__(self):
        tpl_args = {'name': self.name, 'value': _incarify(self.value), 'comment': self.comment, 'comment_sep': self.comment_separator}
        return self._STR_TPL.format(**tpl_args)


class IncarIo(object):
    """
    Parse and write VASP INCAR files.

    Improvements on pymatgen 4.5.3 Incar object:
        * can parse multiple parameters on a line, ';' separated
        * case insensitive
        * parses floats as floats (example: ENCUT truncated to int in pymatgen 4.5.3)

    Limitations:
        * lists of numbers are always parsed as such, in other words a comment after a numeric
            value may not start with a number.

    Writes INCAR files which can be parsed by pymatgen.

    Usage::
        # reuse an INCAR file in an AiiDA calculation
        calc = CalculationFactory('vasp.vasp')
        calc.use_parameters(
           IncarIo(file_path='<path/to/INCAR>').get_param_node()


        # reuse parameter nodes in other tools
        incar_io = IncarIo(parameter_node=calc.inp.parameters)
        ## pymatgen
        incar_io.get_pymatgen()
        ## via string
        str(incar_io)
        ## via file
        incar_io.write('<path/to/file>')

    Parameter name conventions: INCAR reads parameters case insensitively.
    In IncarIO objects parameter names are stored in lower case internally
    and parameter nodes are created with lower case parameter names. This is more
    convenient for typing and ensures query results are consistent.
    When writing to INCAR files, upper case names are used. Pymatgen Incar objects
    do some validation for which upper case names are required.
    """

    def __init__(self, file_path=None, incar_dict=None, parameter_node=None):
        """Populate internal key-value storage from a file or dictionary."""
        self.incar_dict = {}
        if parameter_node:
            self.read_param_node(parameter_node)
        elif file_path:
            self.read_file(file_path)
        elif incar_dict:
            self.incar_dict.update(self.normalize_mapping(incar_dict))

    def read_param_node(self, parameter_node):
        params = parameter_node.get_dict()
        self.incar_dict = self.normalize_mapping(params)

    def read_file(self, file_path):
        """Read an INCAR file into internal storage."""
        with open(file_path) as incar_fo:
            self.read_file_object(incar_fo)

    def read_file_object(self, file_object):
        """Initialize from a file object."""
        self.read_string(file_object.read())

    def read_string(self, input_str):
        """Initialize from a string."""
        self.incar_dict = dict(IncarParamParser.parse_string(input_str))

    @classmethod
    def normalize_mapping(cls, input_mapping):
        normalized_items = [(k.lower(), cls.parse_val(v)) for k, v in input_mapping.items()]
        return dict(normalized_items)

    @classmethod
    def parse_val(cls, input_val):
        if isinstance(input_val, six.string_types):
            return IncarParamParser.value_parser().parseString(input_val)[0]
        return input_val

    def __str__(self):
        items = [IncarItem(*item) for item in sorted(self.incar_dict.items())]
        return '\n'.join((str(incar_item) for incar_item in items))

    def write(self, file_path):
        with open(file_path, 'w') as incar_fo:
            incar_fo.write(str(self))

    def get_dict(self):
        return self.incar_dict

    def get_pymatgen(self):
        return IncarPymatgen.from_dict({k.upper(): v for k, v in self.incar_dict.items()})

    def get_param_node(self):
        """
        Create a ParameterData node containing the incar key/value pairs.

        :kwarg annotate: [True] store the node and add extras to preserve
            order and comments of the INCAR. Implies that the node gets stored in the process!
        """
        node = get_data_node('parameter', dict=self.get_dict())
        return node


class IncarParamParser(object):
    """Parse any INCAR file, including hand written ones."""

    word_chars = pp.alphanums + '-+_:\'"<>~/[]().,'
    item_end = pp.Suppress(pp.MatchFirst(';' | pp.lineEnd))
    equals = pp.Suppress('=')

    @classmethod
    def value_parser(cls):
        """
        Parse INCAR values into an appropriate python type.

        Examples::

            MAGMOM = 1 -1 -> [1, -1]
                     ^^^^ (numeric)

            LORBIT = .False. -> False
                     ^^^^^^^ (bool)

        BNF::

            value               :: numeric | repeated_numeric | boolean | word
            repeated_numeric    :: multiplier+ number
            multiplier          :: integer '*'
            numeric             :: number+
            number              :: integer | decimal
            integer             :: <VASP readable integer literal>
            decimal             :: <VASP readable decimal literal>
            boolean             :: <VASP readable boolean>
        """
        decimal = pp.Regex(r'[+-]?((?:\d+(\.\d*))|(\.\d+))([Ee][+-]?\d+)?').setParseAction(cls.parse_decimal)
        integer = pp.Regex(r'[+-]?(\d+)').setParseAction(cls.parse_int)
        number = decimal | integer
        multiplier = (integer + pp.Literal('*')).setParseAction(cls.parse_int)
        repeated_numeric = (pp.Group(pp.OneOrMore(multiplier)) + number).setParseAction(cls.parse_repeated_num)
        num_value = pp.Group(pp.OneOrMore(number)).setParseAction(cls.parse_num_list)
        false = pp.Regex(r'.f(alse)?.', flags=re.IGNORECASE).setParseAction(lambda t: False)
        true = pp.Regex(r'.t(rue)?.', flags=re.IGNORECASE).setParseAction(lambda t: True)
        bool_value = false | true
        str_value = pp.Word(cls.word_chars)
        return repeated_numeric | num_value | bool_value | str_value

    @classmethod
    def system_parser(cls):
        """
        Parse SYSTEM parameter of INCAR files.

        This is a special case because it can have multi-word values

        Example::

            SYSTEM = This is a system description
            ^----^   ^--------------------------^
             name               value

        BNF::

            system  :: name '=' value
            name    :: 'SYSTEM'
            value   :: (word | ' ' | tab)+
        """
        multiword = pp.Word(cls.word_chars + ' \t')
        multiword_value = multiword + cls.item_end
        system_name = pp.Literal('SYSTEM')
        return system_name('name') + cls.equals + pp.empty + multiword_value('value').setParseAction(lambda t: t[0])

    @classmethod
    def comment_parser(cls):
        """
        Parse comments to INCAR parameter definitions

        Example::

            ENCUT = 280.567 this is an inline comment; LORBIT = .False.
                            ^-----------------------^
                                inline comment

            EMIN = 420 # This is an end line comment; no PARAMETER = assignments after this.
                         ^-----------------------------------------------------------------^
                                    end line comment

        BNF::

            comment         :: (comment_char endline_comment) | inline_comment
            comment_char    :: '#'
            endline_comment :: <everything until end of line>
            inline_comment  :: <anything except ';' or '#' until either ';' or end of line>
        """
        inline_comment = pp.Word(cls.word_chars + ' \t') + cls.item_end
        endline_comment = pp.Suppress('#') + pp.Word(pp.printables + ' \t') + pp.Suppress(pp.lineEnd)
        return endline_comment | inline_comment

    @classmethod
    def get_parser(cls):
        """
        Parse a parameter definition from an INCAR file.

        Realistic examples can be found on
        http://cms.mpi.univie.ac.at/vasp/vasp/INCAR_File.html#incar

        BNF::

            parameter_definition :: system_definition | other_definition
            system_definition    :: system
            other_definition     :: name = value [comment]
        """
        system = cls.system_parser()
        name = pp.Word(pp.alphas, cls.word_chars)
        value = cls.value_parser()
        comment = cls.comment_parser()
        std_param = name('name') + cls.equals + pp.empty + value('value') + pp.Optional(comment)('comment')
        any_param = system | std_param
        param = pp.empty + any_param
        return param

    @classmethod
    def parse_string(cls, input_str):
        all_params = cls.get_parser().scanString(input_str)
        return OrderedDict([(i.name.lower(), i.value) for i, _, _ in all_params])

    @classmethod
    def parse_int(cls, token):
        return int(token[0])

    @classmethod
    def parse_decimal(cls, token):
        return float(token[0])

    @classmethod
    def parse_num_list(cls, token):
        """Parse a token into a list of numbers."""
        num_list = token.asList()
        if len(num_list[0]) == 1:
            return num_list[0]
        return num_list

    @classmethod
    def parse_repeated_num(cls, token):
        multipliers, value = token.asList()
        return {'repetitions': functools.reduce(operator.mul, multipliers), 'value': value}


def _incarify(value):
    """Format value of any compatible type into the string format appropriate for INCAR files."""
    result = None
    if isinstance(value, (str, unicode)):
        result = value
    elif isinstance(value, dict):
        result = "{}*{}".format(value['repetitions'], value['value'])
    elif not np.isscalar(value):
        value_array = np.array(value)
        shape = value_array.shape
        dim = len(shape)
        if dim == 1:
            result = ' '.join([_incarify(i) for i in value])
        elif dim == 2:
            result = '\n'.join([_incarify(i) for i in value])
        elif dim > 2:
            raise TypeError('you are trying to input a more ' + 'than 2-dimensional array to VASP.' + 'Not sure what to do...')
    elif isinstance(value, bool):
        result = '.True.' if value else '.False.'
    elif np.isreal(value):
        result = '{}'.format(value)
    return result


class IncarParser(BaseFileParser):
    """
    Parser for VASP INCAR format.

    This is a wrapper for the parsevasp.incar parser.

    The Parsing direction depends on whether the IncarParser is initialised with
    'path = ...' (read from file) or 'data = ...' (read from data).

    """

    PARSABLE_ITEMS = {
        'incar': {
            'inputs': [],
            'parsers': ['INCAR'],
            'nodeName': '',
            'prerequisites': []
        },
    }

    def __init__(self, *args, **kwargs):
        super(IncarParser, self).__init__(*args, **kwargs)
        self.init_with_kwargs(**kwargs)

    def _init_with_data(self, data):
        """Initialise with a given ParameterData object."""
        self._data_obj = data
        self._parsable_items = self.__class__.PARSABLE_ITEMS
        self._parsed_data = {}

    @property
    def _parsed_object(self):
        """
        Return an instance of parsevasp.incar.Incar.

        Corresponds to the stored data in inputs.parameters.incar.

        """

        incar_dict = self._data_obj.get_dict()

        try:
            return IncarParsevasp(incar_dict=incar_dict, logger=self._logger)
        except SystemExit:
            return None

    def _parse_file(self, inputs):
        """Create a DB Node from an INCAR file."""

        result = inputs
        result = {}

        if isinstance(self._data_obj, get_data_class('parameter')):
            return {'incar': self._data_obj}

        try:
            incar = IncarParsevasp(file_path=self._data_obj.path)
        except SystemExit:
            self._logger.warning("Parsevasp exitited abnormally. Returning None.")
            return {'incar': None}

        result = parsevasp_to_aiida(incar)

        return result


def parsevasp_to_aiida(incar):
    """
    Parsevasp to Aiida conversion.

    Generate an Aiida ParameterData that contains the
    entries found in INCAR using parsevasp.

    """

    incar_dict = incar.get_dict()

    result = {}

    result['incar'] = get_data_class('parameter')(dict=incar_dict)

    return result
