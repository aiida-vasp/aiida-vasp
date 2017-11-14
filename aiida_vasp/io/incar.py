"""Provides INCAR file interface and utilities."""
import re
from collections import OrderedDict

import six
from pymatgen.io.vasp import Incar
import numpy as np
import pyparsing as pp

from aiida_vasp.io.pymatgen_aiida.vasprun import get_data_node
from aiida_vasp.io.parser import KeyValueParser


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
        elif kwargs.has_key('incar_item'):
            self.set_copy(kwargs['incar_item'])
        else:
            self.value = None
            self._name = None
            self._comment = None
            self.set_values(*args, **kwargs)

    @classmethod
    def from_string(cls, input_string):
        result =  cls.item_parser.parseString(input_string)
        value = Incar.proc_val(result.key, result.value)
        return IncarItem(result.key, value, result.comment)

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
        return {
            'name': self.name,
            'value': self.value,
            'comment': self.comment
        }

    def __str__(self):
        tpl_args = {
            'name': self.name,
            'value': _incarify(self.value),
            'comment': self.comment,
            'comment_sep': self.comment_separator
        }
        return self._STR_TPL.format(**tpl_args)


class IncarIo(object):
    """
    Parse and write VASP INCAR files.

    Improvements on pymatgen Incar object:
        * can parse multiple parameters on a line, ';' separated
        * case insensitive

    Writes INCAR files which can be completely parsed by pymatgen.

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
        else:
            if file_path:
                self.read_file(file_path)
            if incar_dict:
                self.update(incar_dict)

    def read_param_node(self, parameter_node):
        params = parameter_node.get_dict()
        self.incar_dict = self.normalize_mapping(params)

    def read_file(self, filename):
        """Read an INCAR file into internal storage."""
        with open(filename) as incar_fo:
            self.read_file_object(incar_fo)

    def read_file_object(self, file_object):
        """Initialize from a file object."""
        self.read_string(file_object.read())

    def read_string(self, input_str):
        """Initialize from a string."""
        self.incar_dict = dict(IncarParamParser.parse_string(input_str))

    @classmethod
    def normalize_mapping(cls, input_mapping):
        normalized_items = [(k.lower(), cls._pmg_proc_val(k, v)) for k, v in input_mapping.items()]
        return dict(normalized_items)

    @classmethod
    def _pmg_proc_val(self, key, input_val):
        if isinstance(input_val, six.string_types):
            return Incar.proc_val(key, input_val)
        return input_val

    def __str__(self):
        items = [IncarItem(*item) for item in sorted(self.incar_dict.items())]
        return '\n'.join((str(incar_item) for incar_item in items))

    def write(self, filename):
        with open(filename, 'w') as incar_fo:
            incar_fo.write(str(self))

    def get_dict(self):
        return self.incar_dict

    def get_pymatgen(self):
        return Incar.from_dict({k.upper(): v for k, v in self.incar_dict.items()})

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
    word_chars = pp.alphanums + '-+_:\'"<>~\/[]().,'
    system_name = pp.Literal('SYSTEM')
    item_end = pp.SkipTo(pp.Suppress(';') | pp.lineEnd)
    equals = pp.Suppress('=')
    system = system_name('name') + equals + item_end('value')
    name = pp.Word(pp.alphas, word_chars)
    number = pp.Regex(r'[+-]?((?:\d+(\.\d*)?)|(.\d+))([Ee][+-]?\d+)?')
    num_value = pp.Group(pp.OneOrMore(number)).setParseAction(lambda t: ' '.join(t[0]))
    str_value = pp.Word(word_chars) + pp.empty
    value = num_value | str_value
    comment = pp.SkipTo(pp.Literal('#') | pp.Literal(';') | pp.lineEnd)
    std_param = name('name') + equals + pp.empty + value('value') + comment('comment')
    any_param = system | std_param
    param = pp.empty + any_param

    @classmethod
    def parse_string(cls, input_str):
        all_params = cls.param.searchString(input_str)
        return OrderedDict([(i.name.lower(), Incar.proc_val(i.name, i.value)) for i in all_params])


def _incarify(value):
    """Format value of any compatible type into the string forat appropriate for INCAR files."""
    result = None
    if isinstance(value, (str, unicode)):
        result = value
    elif not np.isscalar(value):
        value_array = np.array(value)
        shape = value_array.shape
        dim = len(shape)
        if dim == 1:
            result = ' '.join([_incarify(i) for i in value])
        elif dim == 2:
            result = '\n'.join([_incarify(i) for i in value])
        elif dim > 2:
            raise TypeError('you are trying to input a more ' +
                            'than 2-dimensional array to VASP.' +
                            'Not sure what to do...')
    elif isinstance(value, bool):
        result = '.True.' if value else '.False.'
    elif np.isreal(value):
        result = '{}'.format(value)
    return result


def _incar_item(key, value, units=None, comment=None):
    key = key.upper()
    value = _incarify(value)
    units = ' ' + units if units else ''
    comment = ' ' + comment if comment else ''

    return _incar_item.tpl.format(
        key=key, value=value, units=units, comment=comment)


_incar_item.tpl = '{key} = {value}{units}{comment}'


def dict_to_incar(incar_dict, extended=False):
    incar_content = ''
    for key, val in sorted(incar_dict.iteritems(), key=lambda t: t):
        if not extended:
            value = val
            units = None
            comment = None
        else:
            value = val['value']
            units = val.get('units', '')
            comment = val.get('comment', '')
        incar_content += _incar_item(key, value, units, comment) + '\n'
    return incar_content
