"""Utilities for working with aiida in general"""
from functools import wraps
import numpy as np
from packaging import version

from aiida.orm import User

BASIC_DATA_TYPES = ['bool', 'float', 'int', 'list', 'str', 'dict']

def load_dbenv_if_not_loaded(**kwargs):
    """Load dbenv if necessary, run spinner meanwhile to show command hasn't crashed."""
    from aiida.backends.utils import load_dbenv, is_dbenv_loaded
    if not is_dbenv_loaded():
        load_dbenv(**kwargs)


def dbenv(function):
    """A function decorator that loads the dbenv if necessary before running the function."""

    @wraps(function)
    def decorated_function(*args, **kwargs):
        """Load dbenv if not yet loaded, then run the original function."""
        load_dbenv_if_not_loaded()
        return function(*args, **kwargs)

    return decorated_function


def subpath(*args):
    from os.path import dirname, realpath, join
    return realpath(join(dirname(__file__), *args))

def get_data_node(data_type, *args, **kwargs):
    return get_data_class(data_type)(*args, **kwargs)


def get_data_class(data_type):
    """
    Provide access to the orm.data classes with deferred dbenv loading.

    """
    from aiida.plugins import DataFactory
    from aiida.common.exceptions import MissingPluginError

    data_cls = None
    if data_type not in BASIC_DATA_TYPES:
        raise KeyError('Please supply `bool`, `float`, `int`, `list`, `str` or `dict`.')
    try:
        data_cls = DataFactory(data_type)
    except MissingPluginError as err:
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


@dbenv
def create_authinfo(computer, store=False):
    """
    Allow the current user to use the given computer.

    """
    from aiida.orm import AuthInfo
    authinfo = AuthInfo(computer=computer, user=get_current_user())
    if store:
        authinfo.store()
    return authinfo


@dbenv
def cmp_get_authinfo(computer):
    """Get an existing authinfo or None for the given computer and current user."""
    if hasattr(computer, 'get_authinfo'):
        return computer.get_authinfo(get_current_user())
    else:
        from aiida.backends.settings import BACKEND
        from aiida.backends.profile import BACKEND_SQLA, BACKEND_DJANGO

        if BACKEND == BACKEND_DJANGO:
            from aiida.backends.djsite.db.models import DbAuthInfo
            return DbAuthInfo.objects.get(dbcomputer=computer.dbcomputer, aiidauser=get_current_user())  # pylint: disable=no-member
        elif BACKEND == BACKEND_SQLA:
            from aiida.backends.sqlalchemy.models.authinfo import DbAuthInfo
            from aiida.backends.sqlalchemy import get_scoped_session
            session = get_scoped_session()
            return session.query(DbAuthInfo).filter(DbAuthInfo.dbcomputer == computer.dbcomputer).filter(
                DbAuthInfo.aiidauser == get_current_user())
    return None


@dbenv
def cmp_get_transport(computer):
    if hasattr(computer, 'get_transport'):
        return computer.get_transport()
    authinfo = cmp_get_authinfo(computer)
    return authinfo.get_transport()
