"""
Utils for the workchains.

-------------------------
Auxiliary routines that are not part of any of the workchain classes, but needed
to make code more compact in the workchains.
"""
# pylint: disable=import-outside-toplevel
from copy import deepcopy
import numpy as np

from aiida.common.extendeddicts import AttributeDict
from aiida.orm import Dict
from aiida.engine.processes.exit_code import ExitCode
from aiida.plugins import DataFactory

from aiida_vasp.utils.extended_dicts import delete_keys_from_dict


def prepare_process_inputs(inputs, namespaces=None, exclude_parameters=None):
    """
    Prepare the inputs dictionary for a calculation.

    Any remaining bare dictionaries in the inputs dictionary will be wrapped in a Dict data node
    except for the 'options', 'metadata', 'potential' and any key specified in the parameters
    namespaces. They all should remain a standard dictionary. Another exception are dictionaries
    whose keys are not strings but for example tuples.
    """
    from past.builtins import basestring
    prepared_inputs = AttributeDict()

    if namespaces is None:
        namespaces = []

    if exclude_parameters is None:
        exclude_parameters = []

    no_dict = ['options', 'metadata', 'potential', 'parameters']
    no_dict = no_dict + namespaces
    # Copy and convert dict
    for key, val in inputs.items():
        if (key not in no_dict and isinstance(val, dict) and all([isinstance(k, (basestring)) for k in val.keys()])):
            prepared_inputs[key] = Dict(dict=val)
        else:
            prepared_inputs[key] = val

    try:
        # Remove excluded entries for parameters
        parameters = prepared_inputs.parameters
        if exclude_parameters:
            # First make sure we have a proper copy so that any removal does not havoc elements in the dictionary
            if isinstance(parameters, DataFactory('dict')):
                parameters = prepared_inputs.parameters.clone()
                # Unpack in case the parameters is a Dict data structure
                parameters = parameters.get_dict()
            else:
                parameters = deepcopy(prepared_inputs.parameters)
            delete_keys_from_dict(parameters, exclude_parameters)
        if not isinstance(parameters, DataFactory('dict')):
            # Convert parameters to Dict
            parameters = DataFactory('dict')(dict=parameters)
        prepared_inputs.parameters = parameters
    except AttributeError:
        # In case parameters is not present at all
        pass

    return prepared_inputs


def compare_structures(structure_a, structure_b):
    """Compare two StructureData objects A, B and return a delta (A - B) of the relevant properties."""

    delta = AttributeDict()
    delta.absolute = AttributeDict()
    delta.relative = AttributeDict()
    volume_a = structure_a.get_cell_volume()
    volume_b = structure_b.get_cell_volume()
    delta.absolute.volume = np.absolute(volume_a - volume_b)
    delta.relative.volume = np.absolute(volume_a - volume_b) / volume_a

    pos_a = np.array([site.position for site in structure_a.sites])
    pos_b = np.array([site.position for site in structure_b.sites])
    delta.absolute.pos = pos_a - pos_b

    site_vectors = [delta.absolute.pos[i, :] for i in range(delta.absolute.pos.shape[0])]
    a_lengths = np.linalg.norm(pos_a, axis=1)
    delta.absolute.pos_lengths = np.array([np.linalg.norm(vector) for vector in site_vectors])
    delta.relative.pos_lengths = np.array([np.linalg.norm(vector) for vector in site_vectors]) / a_lengths

    cell_lengths_a = np.array(structure_a.cell_lengths)
    delta.absolute.cell_lengths = np.absolute(cell_lengths_a - np.array(structure_b.cell_lengths))
    delta.relative.cell_lengths = np.absolute(cell_lengths_a - np.array(structure_b.cell_lengths)) / cell_lengths_a

    cell_angles_a = np.array(structure_a.cell_angles)
    delta.absolute.cell_angles = np.absolute(cell_angles_a - np.array(structure_b.cell_angles))
    delta.relative.cell_angles = np.absolute(cell_angles_a - np.array(structure_b.cell_angles)) / cell_angles_a

    return delta


def fetch_k_grid(rec_cell, k_spacing):
    """
    Suggest a sensible k-point sampling based on a supplied spacing.

    :param rec_cell: A two dimensional ndarray of floats defining the reciprocal lattice with each vector as row elements.
    :param k_spacing: The k-point spacing.

    :return kgrid: The k-point grid given the supplied `rec_cell` and `kstep`

    This is usable for instance when performing
    plane wave cutoff convergence tests without a base k-point grid.

    """
    rec_cell_lenghts = np.linalg.norm(rec_cell, axis=1)
    kgrid = np.ceil(rec_cell_lenghts / np.float(k_spacing))

    return kgrid.astype('int').tolist()


def compose_exit_code(status, message):
    """Compose an ExitCode instance based on a status and message."""
    exit_code = ExitCode(status=status, message=message)
    return exit_code
