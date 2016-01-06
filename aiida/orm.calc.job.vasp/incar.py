import numpy as np


def _incarify(value):
    result = None
    if isinstance(value, str):
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
        result = value and '.True.' or '.False.'
    elif np.isreal(value):
        result = '{}'.format(value)
    return result


def _incar_item(key, value):
    return _incar_item.tpl.format(key=key.upper(), value=_incarify(value))
_incar_item.tpl = '{key} = {value}'


def dict_to_incar(incar_dict):
    incar_content = ''
    for k, v in incar_dict.iteritems():
        incar_content += _incar_item(k, v) + '\n'
    return incar_content


def par_to_incar(incar_pardat):
    incar_dict = incar_pardat.get_dict()
    return dict_to_incar(incar_dict)
