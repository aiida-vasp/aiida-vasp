"""Provides INCAR file interface and utilities."""
import re
from collections import Mapping, Sequence, OrderedDict

import numpy as np
from six import string_types

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
        self._comment = str(comment).lstrip('#')

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name_input):
        self._name = name.upper()

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

    def __str__(self, rhs):
        return self._STR_TPL.format(**self.as_dict)




class IncarIo(KeyValueParser):
    """Parse and write VASP INCAR files."""
    bool_true = re.compile(r'.true.', re.IGNORECASE)
    bool_false = re.compile(r'.false.', re.IGNORECASE)

    def __init__(self, file_path=None, incar_dict=None):
        """Populate internal key-value storage from a file or dictionary."""
        self.incar_dict = OrderedDict()
        if file_path:
            self.read_incar_file(file_path)
        if incar_dict:
            self.update(incar_dict)

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

        items = re.findall(cls.assignment, content)):
        normalized_items = [(k.lower(), IncarItem(name=k, **cls.clean_value(v))) for k, v in items]
        return OrderedDict(normalized_items)

    @classmethod
    def normalize_mapping(cls, input_mapping):
        normalized_items = [(k.lower(), IncarItem(name=k, **cls.clean_value(v))) for k, v in input_mapping]
        return OrderedDict(normalized_items)

    @classmethod
    def clean_value(cls, str_value):
        cleaned_value = None
        converters = cls.get_converter_iter()
        while not cleaned_value:
            cleaned_value = self.try_convert(converters.next())
        return cleaned_value

    @classmethod
    def get_converter_iter(cls):
        converter_order = [cls.bool, cls.int, cls.float, cls.string]
        return (i for i in converter_order)

    def try_convert(cls, str_value, converter):
        try:
            cleaned_value = converter(str_value)
        except ValueError:
            cleaned_value = {}

        if cleaned_value.get('value', None) is None:
            return None
        return cleaned_value

    def update(self, incar_mapping):
        input_dict = self.normalize_mapping(incar_mapping)
        self.incar_dict.update(input_dict)

    def __str__(self):
        return dict_to_incar(self.incar_dict, extended=True)

    def store(self, filename):
        with open(filename, 'w') as incar_fo:
            incar_fo.write(str(self))

    def get_dict(self):
        return {v.name: v.value for v in self.incar_dict.values()}

    def make_param_node(self):
        return get_data_node('param', dict=self.get_dict())



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

    return _incar_item.tpl.format(key=key, value=value, units=units, comment=comment)


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
