"""
Utils for AiiDA.

----------------
Utilities for making working against AiiDA a bit easier. Mostly here due to
historical reasons when AiiDA was rapidly developed. In the future most routines
that have now standardized in AiiDA will be removed.
"""

import warnings
from functools import wraps

import numpy as np
from aiida import orm
from aiida.orm import User, load_node
from packaging import version

BASIC_DATA_TYPES = ['core.bool', 'core.float', 'core.int', 'core.list', 'core.str', 'core.dict']


def querybuild(cls, **kwargs):
    """
    Instantiates and returns a QueryBuilder instance.

    The QueryBuilder's path has one vertice so far, namely this class.
    Additional parameters (e.g. filters or a label),
    can be passes as keyword arguments.

    :param label: Label to give
    :param filters: filters to apply
    :param project: projections
    :returns: a QueryBuilder instance.
    """

    from aiida.orm import QueryBuilder

    query_builder = QueryBuilder()
    filters = kwargs.pop('filters', {})
    query_builder.append(cls, filters=filters, **kwargs)

    return query_builder


def get_data_class(data_type):
    """Provide access to the orm.data classes with deferred dbenv loading."""
    from aiida.common.exceptions import MissingEntryPointError
    from aiida.plugins import DataFactory

    data_cls = None
    try:
        data_cls = DataFactory(data_type)
    except MissingEntryPointError as err:
        raise err
    return data_cls


def get_current_user():
    """Get current user."""
    current_user = User.collection.get_default()
    return current_user


def copy_parameter(old_parameter):
    """Assemble a new Dict."""
    return orm.Dict(dict=old_parameter.get_dict())


def displaced_structure(structure, displacement, entry):
    disp_structure = structure.clone()
    displace_position(disp_structure, displacement, entry)
    return disp_structure


def compressed_structure(structure, volume_change):
    comp_structure = structure.clone()
    compress_cell(comp_structure, volume_change)
    return comp_structure


def displace_position(structure, displacement, entry):
    """Displace a position in the StructureData."""
    sites = structure.sites
    positions = []
    for site in sites:
        positions.append(site.position)
    new_position = np.asarray(positions[entry - 1]) + displacement
    new_position = new_position.tolist()
    positions[entry - 1] = tuple(new_position)
    structure.reset_sites_positions(positions)


def compress_cell(structure, volume_change):
    """Apply compression or tensile forces to the unit cell."""
    cell = structure.cell
    new_cell = np.array(cell) * volume_change
    structure.reset_cell(new_cell.tolist())


def aiida_version():
    from aiida import __version__ as aiida_version_

    return version.parse(aiida_version_)


def cmp_version(string):
    return version.parse(string)


def cmp_load_verdi_data():
    """Load the verdi data click command group for any version since 0.11."""
    verdi_data = None
    import_errors = []

    try:
        from aiida.cmdline.commands import data_cmd as verdi_data
    except ImportError as err:
        import_errors.append(err)

    if not verdi_data:
        try:
            from aiida.cmdline.commands import verdi_data
        except ImportError as err:
            import_errors.append(err)

    if not verdi_data:
        try:
            from aiida.cmdline.commands.cmd_data import verdi_data
        except ImportError as err:
            import_errors.append(err)

    if not verdi_data:
        err_messages = '\n'.join([f' * {err}' for err in import_errors])
        raise ImportError('The verdi data base command group could not be found:\n' + err_messages)

    return verdi_data


def create_authinfo(computer, store=False):
    """Allow the current user to use the given computer."""
    from aiida.orm import AuthInfo

    authinfo = AuthInfo(computer=computer, user=get_current_user())
    if store:
        authinfo.store()
    return authinfo


def cmp_get_authinfo(computer):
    """Get an existing authinfo or None for the given computer and current user."""
    return computer.get_authinfo(get_current_user())


def cmp_get_transport(computer):
    if hasattr(computer, 'get_transport'):
        return computer.get_transport()
    authinfo = cmp_get_authinfo(computer)
    return authinfo.get_transport()


def ensure_node_first_arg(func):
    """Decorator to load a node if it is passed as a string."""

    @wraps(func)
    def wrapper(*args, **kwargs):
        """Make sure the first node is a Node instance."""
        if len(args) > 0:
            node = args[0]
            if not isinstance(node, orm.Node):
                node = load_node(node)
        args = list(args)
        args[0] = node
        return func(*args, **kwargs)

    return wrapper


def ensure_node_kwargs(func):
    """Decorator to load a node if it is passed as a key word argument ends with 'node'."""

    @wraps(func)
    def wrapper(node, *args, **kwargs):
        """Make sure the key world arguments ends with '_node' node is a Node instance."""
        new_kwargs = dict(kwargs)
        for name, value in kwargs.items():
            if name.endswith('node'):
                if not isinstance(kwargs[name], orm.Node):
                    new_kwargs[name] = load_node(value)
        return func(node, *args, **new_kwargs)

    return wrapper


def convert_dict_case(dict_in: dict, recursive=True, warn=False, lower=True, raise_convert=False):
    """
    Recursively convert the keys of a dictionary to lower or upper cases, returns a new dictionary.

    :param dict_in: The input dictionary whose keys need to be converted.
    :param recursive: If True, the function will recursively convert keys in nested dictionaries.
    :param warn: If True, the function will print a warning if a key is converted.
    :param lower: If True, convert keys to lowercase; otherwise, convert to uppercase.
    :param raise_convert: If True, raise an error if a key is converted.
    :return: A new dictionary with keys converted to the specified case.
    """

    converted_dict = {}
    for key, value in dict_in.items():
        new_key = key.lower() if lower else key.upper()
        if new_key != key:
            expected = 'upper' if lower is False else 'lower'
            if warn:
                expected = 'upper' if lower is False else 'lower'
                warnings.warn(f"Key '{key}' converted to '{new_key}' - please use {expected} case keys")
            if raise_convert:
                raise ValueError(f"Key '{key}' converted to '{new_key}' - please use {expected} case keys")

        if recursive and isinstance(value, dict):
            converted_dict[new_key] = convert_dict_case(value, recursive, warn, lower, raise_convert)
        else:
            converted_dict[new_key] = value
    return converted_dict
