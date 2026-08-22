"""
Common functions and constants
"""

# Name of the override name space
# This is the namespace where raw VASP INCAR tags should reside for VaspWorkChain
import warnings
from functools import wraps
from typing import Any, Callable

from aiida import orm
from aiida.common.exceptions import InputValidationError
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.assistant.parameters import _BASE_NAMESPACES, ParametersMassage
from aiida_vasp.utils.aiida_utils import convert_dict_case

OVERRIDE_NAMESPACE = 'incar'

# pylint:disable=raise-missing-from


def aiida_to_python(entity: Any) -> Any:
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


def plain_python_args(func: Callable[..., Any]) -> Callable[..., Any]:
    """Ensure that the first argument is a plain dictionary"""

    @wraps(func)
    def wrapped(*args: Any, **kwargs: Any) -> Any:
        new_args = list(args)
        new_args[0] = aiida_to_python(args[0])
        return func(*new_args, **kwargs)

    return wrapped


def parameters_validator(node: orm.Dict | None, port: Any = None) -> None:
    """
    Validate the parameters input by passing it through the massager
    """
    _ = port
    if not node:
        return

    pdict = node.get_dict()

    try:
        convert_dict_case(pdict, lower=True, raise_convert=node.is_stored)
    except ValueError as error:
        raise InputValidationError(
            f'Case inconsistency found in the parameters dictionary please use lower case keys: {error}'
        )

    if not node.is_stored:
        node.set_dict(convert_dict_case(pdict, lower=True, warn=True))

    if OVERRIDE_NAMESPACE not in pdict:
        raise InputValidationError(f'Would expect some incar tags supplied under {OVERRIDE_NAMESPACE} key!')

    accepted_namespaces = _BASE_NAMESPACES + [OVERRIDE_NAMESPACE]
    new_dict = {key: value for key, value in pdict.items() if key in accepted_namespaces}
    try:
        ParametersMassage(new_dict)
    except Exception as error:
        raise InputValidationError(f'Cannot validate the input parameters - error from massager: {error}')


def warn_deprecated_options(node: orm.Dict | None, port: Any = None) -> None:
    """
    Validate the parameters input by passing it through the massager
    """
    _ = port
    _ = node
    warnings.warn('The use of `options` port is deprecated, please use `calc.metadata.options` instead.')


@plain_python_args
def site_magnetization_to_magmom(site_dict: dict[str, Any]) -> list[Any]:
    """
    Convert site mangetization to MAGMOM used for restart.

    For collinear calculations (ISPIN = 2) the function returns a list of
    scalar magnetic moments per site. For non-collinear calculations
    (LNONCOLLINEAR = .TRUE.) the function returns a list of 3-tuples
    ``(mx, my, mz)`` per site.

    NOTE: to be replaced by stock function in aiida_vasp.utils.workchains
    """
    if 'site_magnetization' in site_dict:
        site_dict = site_dict['site_magnetization']

    site_dict = site_dict['sphere']

    # Detect whether we are dealing with a non-collinear case by checking
    # which directions are populated.
    available_axes = []
    for symbol in 'xyz':
        if site_dict.get(symbol) and site_dict.get(symbol, {}).get('site_moment'):
            available_axes.append(symbol)

    if not available_axes:
        raise ValueError('No valid site-projected magnetization available')

    # Choose the axes used to build the returned magmom. For collinear only
    # one axis is populated; for noncollinear we use all available axes.
    if set(available_axes) >= {'x', 'y', 'z'} or len(available_axes) > 1:
        axes = [axis for axis in 'xyz' if axis in available_axes]
    else:
        axes = available_axes

    # Ensure each axis is sorted by site index
    sorted_axes = {}
    for axis in axes:
        tmp = list(site_dict[axis]['site_moment'].items())
        tmp.sort(key=lambda x: int(x[0]))
        sorted_axes[axis] = [entry[1]['tot'] for entry in tmp]

    n_sites = len(sorted_axes[axes[0]])
    if len(axes) == 1:
        return list(sorted_axes[axes[0]])
    # Non-collinear: return one 3-tuple per site (zeros are inserted for
    # missing axes)
    magmom: list[Any] = []
    for i in range(n_sites):
        components = [float(sorted_axes[axis][i]) if axis in sorted_axes else 0.0 for axis in 'xyz']
        magmom.append(tuple(components))
    return magmom


def nested_update(dict_in: dict[str, Any], update_dict: dict[str, Any], extend_list: bool = False) -> dict[str, Any]:
    """Update the dictionary - combine nested sub-dictionary with update as well"""
    warnings.warn('nested_update is deprecated, use updated_nested_dict', DeprecationWarning)
    for key, value in update_dict.items():
        if key in dict_in and isinstance(value, (dict, AttributeDict)):
            nested_update(dict_in[key], value, extend_list=extend_list)
        elif key in dict_in and isinstance(value, list) and extend_list:
            dict_in[key].extend(value)
        else:
            dict_in[key] = value
    return dict_in


def nested_update_dict_node(dict_node: orm.Dict, update_dict: dict[str, Any], extend_list: bool = False) -> orm.Dict:
    """Utility to update a Dict node in a nested way"""
    warnings.warn('nested_update_dict_node is deprecated, use updated_nested_dict_node', DeprecationWarning)
    pydict = dict_node.get_dict()
    nested_update(pydict, update_dict, extend_list=extend_list)
    if pydict == dict_node.get_dict():
        return dict_node
    return orm.Dict(dict=pydict)
