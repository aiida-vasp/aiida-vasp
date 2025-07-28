"""
Utils for AiiDA.

Utilities for making working against AiiDA a bit easier. Mostly here due to
historical reasons when AiiDA was rapidly developed. In the future most routines
that have now standardized in AiiDA will be removed.
"""

# ruff: noqa: PLC0415
from __future__ import annotations

import warnings
from functools import wraps
from typing import Any, Callable

import numpy as np
from aiida import __version__ as aiida_version_
from aiida import orm
from aiida.common.exceptions import MissingEntryPointError
from aiida.orm import AuthInfo, QueryBuilder, User, load_node
from aiida.plugins import DataFactory
from packaging import version

BASIC_DATA_TYPES: list[str] = ['core.bool', 'core.float', 'core.int', 'core.list', 'core.str', 'core.dict']


def querybuild(cls: type, **kwargs: Any) -> QueryBuilder:
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

    query_builder = QueryBuilder()
    filters = kwargs.pop('filters', {})
    query_builder.append(cls, filters=filters, **kwargs)

    return query_builder


def get_data_class(data_type: str) -> type:
    """Provide access to the orm.data classes with deferred dbenv loading."""

    data_cls = None
    try:
        data_cls = DataFactory(data_type)
    except MissingEntryPointError as err:
        raise err
    return data_cls


def get_current_user() -> User:
    """Get current user."""
    current_user = User.collection.get_default()
    return current_user


def copy_parameter(old_parameter: orm.Dict) -> orm.Dict:
    """Assemble a new Dict."""
    return orm.Dict(dict=old_parameter.get_dict())


def displaced_structure(structure: orm.StructureData, displacement: np.ndarray, entry: int) -> orm.StructureData:
    disp_structure = structure.clone()
    displace_position(disp_structure, displacement, entry)
    return disp_structure


def compressed_structure(structure: orm.StructureData, volume_change: float) -> orm.StructureData:
    comp_structure = structure.clone()
    compress_cell(comp_structure, volume_change)
    return comp_structure


def displace_position(structure: orm.StructureData, displacement: np.ndarray, entry: int) -> None:
    """Displace a position in the StructureData."""
    sites = structure.sites
    positions = []
    for site in sites:
        positions.append(site.position)
    new_position = np.asarray(positions[entry - 1]) + displacement
    new_position = new_position.tolist()
    positions[entry - 1] = tuple(new_position)
    structure.reset_sites_positions(positions)


def compress_cell(structure: orm.StructureData, volume_change: float) -> None:
    """Apply compression or tensile forces to the unit cell."""
    cell = structure.cell
    new_cell = np.array(cell) * volume_change
    structure.reset_cell(new_cell.tolist())


def aiida_version() -> version.Version:
    return version.parse(aiida_version_)


def cmp_version(string: str) -> version.Version:
    return version.parse(string)


def cmp_load_verdi_data() -> Any:
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


def create_authinfo(computer: orm.Computer, store: bool = False) -> AuthInfo:
    """Allow the current user to use the given computer."""

    authinfo = AuthInfo(computer=computer, user=get_current_user())
    if store:
        authinfo.store()
    return authinfo


def cmp_get_authinfo(computer: orm.Computer) -> AuthInfo | None:
    """Get an existing authinfo or None for the given computer and current user."""
    return computer.get_authinfo(get_current_user())


def cmp_get_transport(computer: orm.Computer) -> Any:
    if hasattr(computer, 'get_transport'):
        return computer.get_transport()
    authinfo = cmp_get_authinfo(computer)
    return authinfo.get_transport()


def ensure_node_first_arg(func: Callable[..., Any]) -> Callable[..., Any]:
    """Decorator to load a node if it is passed as a string."""

    @wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        """Make sure the first node is a Node instance."""
        if len(args) > 0:
            node = args[0]
            if not isinstance(node, orm.Node):
                node = load_node(node)
        args = list(args)
        args[0] = node
        return func(*args, **kwargs)

    return wrapper


def ensure_node_kwargs(func: Callable[..., Any]) -> Callable[..., Any]:
    """Decorator to load a node if it is passed as a key word argument ends with 'node'."""

    @wraps(func)
    def wrapper(node: Any, *args: Any, **kwargs: Any) -> Any:
        """Make sure the key world arguments ends with '_node' node is a Node instance."""
        new_kwargs = dict(kwargs)
        for name, value in kwargs.items():
            if name.endswith('node'):
                if not isinstance(value, orm.Node):
                    new_kwargs[name] = load_node(value)
        return func(node, *args, **new_kwargs)

    return wrapper


def convert_dict_case(
    dict_in: dict[str, Any],
    recursive: bool = True,
    warn: bool = False,
    lower: bool = True,
    raise_convert: bool = False,
) -> dict[str, Any]:
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
