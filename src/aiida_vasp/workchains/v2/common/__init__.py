"""
Common functions and constants
"""
# Name of the override name space
# This is the namespace where raw VASP INCAR tags should reside for VaspWorkChain
from functools import wraps

from aiida import orm
from aiida.common.exceptions import InputValidationError
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.assistant.parameters import _BASE_NAMESPACES, ParametersMassage

OVERRIDE_NAMESPACE = 'incar'

#pylint:disable=raise-missing-from


def aiida_to_python(entity):
    """
    Convert AiiDA entity to plain python objects
    """
    if not isinstance(entity, orm.Data):
        return entity
    if isinstance(entity, orm.Dict):
        return entity.get_dict()
    if isinstance(entity, orm.List):
        return entity.get_list()
    if isinstance(entity, (orm.Float, orm.Str, orm.Int)):
        return entity.value
    raise ValueError(f'{entity} cannot be converted to plain python object')


def plain_python_args(func):
    """Ensure that the first argument is a plain dictionary"""

    @wraps(func)
    def wrapped(*args, **kwargs):
        new_args = list(args)
        new_args[0] = aiida_to_python(args[0])
        return func(*new_args, **kwargs)

    return wrapped


def parameters_validator(node, port=None):
    """
    Validate the parameters input by passing it through the massager
    """
    _ = port
    if not node:
        return

    pdict = node.get_dict()
    if OVERRIDE_NAMESPACE not in pdict:
        raise InputValidationError(f'Would expect some incar tags supplied under {OVERRIDE_NAMESPACE} key!')

    accepted_namespaces = _BASE_NAMESPACES + [OVERRIDE_NAMESPACE]
    new_dict = {key: value for key, value in pdict.items() if key in accepted_namespaces}
    try:
        ParametersMassage(new_dict)
    except Exception as error:
        raise InputValidationError(f'Cannot validate the input parameters - error from massager: {error}')


@plain_python_args
def site_magnetization_to_magmom(site_dict):
    """
    Convert site mangetization to MAGMOM used for restart
    NOTE: to be replaced by stock function in aiida_vasp.utils.workchains
    """
    if 'site_magnetization' in site_dict:
        site_dict = site_dict['site_magnetization']

    site_dict = site_dict['sphere']
    to_use = None
    for symbol in 'xyz':
        if site_dict.get(symbol) and site_dict.get(symbol, {}).get('site_moment'):
            to_use = symbol
            break
    # No avaliable site magnetization for setting MAGMOM, something is wrong
    if to_use is None:
        raise ValueError('No valid site-projected magnetization avaliable')
    # Ensure sorted list
    tmp = list(site_dict[to_use]['site_moment'].items())
    tmp.sort(key=lambda x: int(x[0]))
    return [entry[1]['tot'] for entry in tmp]


def nested_update(dict_in, update_dict):
    """Update the dictionary - combine nested subdictionary with update as well"""
    for key, value in update_dict.items():
        if key in dict_in and isinstance(value, (dict, AttributeDict)):
            nested_update(dict_in[key], value)
        else:
            dict_in[key] = value
    return dict_in


def nested_update_dict_node(dict_node, update_dict):
    """Utility to update a Dict node in a nested way"""
    pydict = dict_node.get_dict()
    nested_update(pydict, update_dict)
    if pydict == dict_node.get_dict():
        return dict_node
    return orm.Dict(dict=pydict)
