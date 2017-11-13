"""Provides INCAR file interface and utilities."""
import re
from collections import OrderedDict

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


class IncarIo(KeyValueParser):
    """Parse and write VASP INCAR files, preserving order and comments."""
    bool_true = re.compile(r'.true.', re.IGNORECASE)
    bool_false = re.compile(r'.false.', re.IGNORECASE)
    key_order = '__order'
    key_comments = '__comments'

    def __init__(self, file_path=None, incar_dict=None, parameter_node=None):
        """Populate internal key-value storage from a file or dictionary."""
        self.incar_dict = OrderedDict()
        if parameter_node:
            self.read_param_node(parameter_node)
        else:
            if file_path:
                self.read_incar_file(file_path)
            if incar_dict:
                self.update(incar_dict)

    def read_param_node(self, parameter_node):
        params = parameter_node.get_dict()
        order = params.pop(self.key_order, [])
        comments = params.pop(self.key_comments, {})
        self.incar_dict = OrderedDict([(k, None) for k in order])
        self.incar_dict.update({k: IncarItem(k, v, comments.get(k, None)) for k, v in params.items()})

    def add_defaults(self, defaults_mapping):
        for k, v in self.normalize_mapping(defaults_mapping):
            if k not in self.incar_dict:
                self.incar_dict[k] = v

    def read_incar_file(self, filename):
        """Read an INCAR file into internal storage."""
        with open(filename) as incar:
            self.incar_dict = self.parse_incar_fo(incar)

    @classmethod
    def parse_incar_fo(cls, fobj_or_str):
        """Read key/value pairs from a file object into a dictionary."""
        if isinstance(fobj_or_str, str):
            content = fobj_or_str
        else:
            content = fobj_or_str.read()

        items = re.findall(cls.assignment, content)
        normalized_items = [(k.lower(), IncarItem(
            name=k, **cls.clean_value(v))) for k, v in items]
        return OrderedDict(normalized_items)

    @classmethod
    def normalize_mapping(cls, input_mapping):
        normalized_items = [(k.lower(), IncarItem(
            name=k, **cls.clean_value(v))) for k, v in input_mapping.items()]
        return OrderedDict(normalized_items)

    def update(self, incar_mapping):
        input_dict = self.normalize_mapping(incar_mapping)
        self.incar_dict.update(input_dict)

    def __str__(self):
        return '\n'.join((str(incar_item)
                          for incar_item in self.incar_dict.values()))

    def store(self, filename):
        with open(filename, 'w') as incar_fo:
            incar_fo.write(str(self))

    def get_dict(self):
        return {k: v.value for k, v in self.incar_dict.items()}

    def get_comments_dict(self):
        return {k: v.comment for k, v in self.incar_dict.items()}

    def make_param_node(self, annotate=True):
        """
        Create a ParameterData node containing the incar key/value pairs.

        :kwarg annotate: [True] store the node and add extras to preserve
            order and comments of the INCAR. Implies that the node gets stored in the process!
        """

        node = get_data_node('parameter', dict=self.get_dict())
        if not annotate:
            return node
        node.store()
        node.set_extra(self.extra_key_order, self.incar_dict.keys())
        node.set_extra(self.extra_key_comments, self.get_comments_dict())
        return node



def item_parser():
    name_token = pp.Word(pp.alphas, pp.alphanums + '_')
    assignment_token = pp.Suppress('=')
    nums_token = pp.Regex(r'(?P<list>(?P<num>[\-\+]?[\d.,\']+e?[\-\+]?\d*)\s?)*').setParseAction(lambda t: t[0].strip())
    value_token = nums_token | pp.Word(pp.alphanums + '_.,+-')
    item_end_token = pp.SkipTo(pp.Suppress(';') | pp.lineEnd)
    comment_token = pp.Suppress(pp.Optional('#')) + item_end_token('comment')
    item_token = name_token('key') + assignment_token + value_token('value') + comment_token
    description_key_token = pp.Literal('SYSTEM')
    description_token = description_key_token('key') + assignment_token + item_end_token('value')
    item = description_token | item_token
    return item


def item_parser2():
    word_chars = pp.printables.replace(';', '')
    system_name = pp.Literal('SYSTEM')
    item_end = pp.SkipTo(pp.Suppress(';') | pp.lineEnd)
    equals = pp.Suppress('=')
    system = system_name('name') + equals + item_end('value')
    name = pp.Word(pp.alphas, word_chars)
    number = pp.Regex(r'[+-]?((?:\d+(\.\d*)?)|(.\d+))([Ee][+-]?\d+)?')
    num_value = pp.Group(pp.OneOrMore(number)).setParseAction(lambda t: ' '.join(t[0]))
    str_value = pp.Word(word_chars)
    value = num_value | str_value
    comment = pp.SkipTo(pp.Literal('#') | pp.Literal(';') | pp.lineEnd)
    std_param = name('name') + equals + pp.empty + value('value') + comment('comment')
    param = system | std_param
    return param

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
