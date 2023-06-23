"""
Utils for AiiDA.

----------------
Utilities for making working against AiiDA a bit easier. Mostly here due to
historical reasons when AiiDA was rapidly developed. In the future most routines
that have now standardized in AiiDA will be removed.
"""
# pylint: disable=import-outside-toplevel
import numpy as np
from packaging import version

from aiida.orm import User

from aiida.cmdline.utils.decorators import with_dbenv

BASIC_DATA_TYPES = ['bool', 'float', 'int', 'list', 'str', 'dict']


def get_data_node(data_type, *args, **kwargs):
    return get_data_class(data_type)(*args, **kwargs)


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


@with_dbenv()
def get_data_class(data_type):
    """Provide access to the orm.data classes with deferred dbenv loading."""
    from aiida.plugins import DataFactory
    from aiida.common.exceptions import MissingEntryPointError

    data_cls = None
    try:
        data_cls = DataFactory(data_type)
    except MissingEntryPointError as err:
        raise err
    return data_cls


def get_current_user():
    """Get current user."""
    current_user = User.objects.get_default()
    return current_user


def copy_parameter(old_parameter):
    """Assemble a new Dict."""
    return get_data_node('dict', dict=old_parameter.get_dict())


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
        err_messages = '\n'.join([' * {}'.format(err) for err in import_errors])
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
