"""Provides INCAR file interface and utilities."""
import re

import numpy as np

from aiida_vasp.io.parser import KeyValueParser


class IncarIo(KeyValueParser):
    """Parse and write VASP INCAR files."""
    bool_true = re.compile(r'.true.', re.IGNORECASE)
    bool_false = re.compile(r'.false.', re.IGNORECASE)

    def __init__(self, file_path=None, incar_dict=None):
        """Populate internal key-value storage from a file or dictionary."""
        self.incar_dict = {}
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
        kvd = dict(re.findall(cls.assignment, content))
        kvd = {k.lower(): cls.clean_value(v) for k, v in kvd.iteritems()}
        return kvd

    @classmethod
    def clean_value(cls, str_value):
        cleaned_value = {'value': None}
        try:
            cleaned_value = cls.bool(str_value)
        except ValueError as err:
            pass

        if cleaned_value['value'] is None:
            try:
                cleaned_value = cls.int_unit(str_value)
            except ValueError:
                pass

        if cleaned_value['value'] is None:
            try:
                cleaned_value = cls.float_unit(str_value)
            except ValueError:
                pass

        if cleaned_value['value'] is None:
            cleaned_value = cls.string(str_value)

        return cleaned_value

    def update(self, incar_dict):
        # TODO: This is really ugly and the whole storage format should
        # TODO: be changed to make this not so weird
        for key, val in incar_dict.items():
            if key not in self.incar_dict:
                self.incar_dict[key.lower()] = {}
            if not hasattr(val, 'has_key'):
                self.incar_dict[key.lower()]['value'] = val
            else:
                self.incar_dict[key.lower()].update(val)

    def __str__(self):
        return dict_to_incar(self.incar_dict, extended=True)

    def store(self, filename):
        with open(filename, 'w') as incar_fo:
            incar_fo.write(str(self))



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
