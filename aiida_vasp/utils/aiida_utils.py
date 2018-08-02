"""Utilities for working with aiida in general"""
from functools import wraps
import numpy as np

from aiida.common.extendeddicts import AttributeDict


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


def init_input(inputs):
    """
    Assemble the input into a AttributeDict.

    It is necessary to do this if passing the inputs come from a previous
    submit and one is set to pass it along.
    """
    assembled_inputs = AttributeDict()
    for name, value in inputs.items():
        assembled_inputs[name] = value

    return assembled_inputs


def new_structure(old_structure):
    """Assemble a new StructureData."""
    structure_cls = get_data_class('structure')
    old_cell = old_structure.cell
    structure = structure_cls(cell=old_cell)
    old_sites = old_structure.sites
    kinds = old_structure.kinds
    for index, site in enumerate(old_sites):
        structure.append_atom(position=site.position, symbols=kinds[index].symbols, name=kinds[index].name)
    return structure


def new_parameter(old_parameter):
    """Assemble a new ParameterData."""
    parameter_cls = get_data_class('parameter')
    return parameter_cls(dict=old_parameter.get_dict())


def new_kpoints(old_kpoints, structure):
    """Assemble a new KpointsData."""
    kpoints = get_data_class('array.kpoints')
    mesh = old_kpoints.get_kpoints_mesh()
    return kpoints(kpoints=mesh, cell_from_structure=structure)


def displaced_structure(structure, displacement, entry):
    disp_structure = new_structure(structure)
    displace_position(disp_structure, displacement, entry)
    return disp_structure


def compressed_structure(structure, volume_change):
    comp_structure = new_structure(structure)
    compress_cell(comp_structure, volume_change)
    return comp_structure


def displace_position(structure, displacement, entry):
    """Displace a position in the StructureData."""
    sites = structure.sites
    positions = []
    for site in sites:
        positions.append(site.position)
    position = positions[entry]
    position = position + displacement
    structure.reset_sites_positions(positions)


def compress_cell(structure, volume_change):
    """Apply compression or tensile forces to the unit cell."""
    cell = structure.cell
    new_cell = np.array(cell) * volume_change
    structure.reset_cell(new_cell.to_list())
