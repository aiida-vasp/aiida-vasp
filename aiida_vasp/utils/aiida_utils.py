"""Utilities for working with aiida in general"""
from functools import wraps
from packaging import version


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


def get_data_node(data_type, *args, **kwargs):
    return get_data_class(data_type)(*args, **kwargs)


@dbenv
def get_data_class(data_type):
    """
    Provide access to the orm.data classes with deferred dbenv loading.

    compatiblity: also provide access to the orm.data.base memebers, which are loadable through the DataFactory as of 1.0.0-alpha only.
    """
    from aiida.orm import DataFactory
    from aiida.common.exceptions import MissingPluginError
    data_cls = None
    try:
        data_cls = DataFactory(data_type)
    except MissingPluginError as err:
        if data_type in BASIC_DATA_TYPES:
            data_cls = get_basic_data_pre_1_0(data_type)
        else:
            raise err
    return data_cls


BASIC_DATA_TYPES = set(['bool', 'float', 'int', 'list', 'str'])


@dbenv
def get_basic_data_pre_1_0(data_type):
    from aiida.orm.data import base as base_data
    return getattr(base_data, data_type.capitalize())


@dbenv
def backend_obj_users():
    """Test if aiida accesses users through backend object."""
    backend_obj_flag = False
    try:
        from aiida.backends.utils import get_automatic_user  # pylint: disable=unused-variable,no-name-in-module
    except ImportError:
        backend_obj_flag = True
    return backend_obj_flag


@dbenv
def get_current_user():
    """Get current user backwards compatibly with aiida-core <= 0.12.1."""
    current_user = None
    if backend_obj_users():
        from aiida.orm.backend import construct_backend  # pylint: disable=no-name-in-module
        backend = construct_backend()
        current_user = backend.users.get_automatic_user()
    else:
        from aiida.backends.utils import get_automatic_user  # pylint: disable=no-name-in-module
        current_user = get_automatic_user()
    return current_user


def builder_interface(calc_cls):
    """Return the JobProcess or the JobCalculation class, depending on aiida version."""
    if hasattr(calc_cls, 'get_builder'):
        return True
    return False


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

    Deal with backwards compatibility down to aiida 0.11
    """
    from aiida.orm import backend as orm_backend
    authinfo = None
    if hasattr(orm_backend, 'construct_backend'):
        backend = orm_backend.construct_backend()
        authinfo = backend.authinfos.create(computer=computer, user=get_current_user())
        if store:
            authinfo.store()
    else:
        from aiida.backends.settings import BACKEND
        from aiida.backends.profile import BACKEND_SQLA, BACKEND_DJANGO

        if BACKEND == BACKEND_DJANGO:
            from aiida.backends.djsite.db.models import DbAuthInfo
            authinfo = DbAuthInfo(dbcomputer=computer.dbcomputer, aiidauser=get_current_user())
        elif BACKEND == BACKEND_SQLA:
            from aiida.backends.sqlalchemy.models.authinfo import DbAuthInfo
            from aiida.backends.sqlalchemy import get_scoped_session
            _ = get_scoped_session()
            authinfo = DbAuthInfo(dbcomputer=computer.dbcomputer, aiidauser=get_current_user())
        if store:
            authinfo.save()
    return authinfo


@dbenv
def cmp_get_transport(computer, authinfo):
    if hasattr(computer, 'get_transport'):
        return computer.get_transport()
    return authinfo.get_transport()
