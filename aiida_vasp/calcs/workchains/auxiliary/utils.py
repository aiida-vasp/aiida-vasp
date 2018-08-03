"""Auxiliary routines that are not part of any of the workchain classes."""

import numpy as np
from aiida.common.extendeddicts import AttributeDict
from aiida.orm.data.parameter import ParameterData


def prepare_process_inputs(inputs):
    """
    Prepare the inputs dictionary for a calculation.

    Any remaining bare dictionaries in the inputs dictionary will be wrapped in a ParameterData data node
    except for the '_options/options' key which should remain a standard dictionary. Another exception are dictionaries
    whose keys are not strings but for example tuples.
    This is the format used by input groups as in for example the explicit pseudo dictionary where the key is
    a tuple of kind to which the UpfData corresponds.
    """
    prepared_inputs = AttributeDict()

    for key, val in inputs.iteritems():
        if key not in ['_options', 'options'] and isinstance(val, dict) and all([isinstance(k, (basestring)) for k in val.keys()]):
            prepared_inputs[key] = ParameterData(dict=val)
        else:
            prepared_inputs[key] = val

    return prepared_inputs


def finished_ok_compat(calc):
    if hasattr(calc, 'has_finished_ok'):
        return calc.has_finished_ok()
    elif hasattr(calc, 'is_finished_ok'):
        return calc.is_finished_ok
    return calc.finished_ok


def compare_structures(structure_a, structure_b):
    """Compare two StructreData objects A, B and return a delta (A - B) of the relevant properties."""

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

    Parameters
    ----------
    rec_cell : ndarray
        A two dimensional ndarray of floats defining the reciprocal lattice with each
        vector as row elements.
    k_spacing : float
        The k-point spacing.

    Returns
    -------
    kgrid : (3) list of int
        The k-point grid given the supplied `rec_cell` and `kstep`

    Notes
    -----
    This is usable for instance when performing
    plane wave cutoff convergence tests without a base k-point grid.

    """
    rec_cell_lenghts = np.linalg.norm(rec_cell, axis=1)
    kgrid = np.ceil(rec_cell_lenghts / np.float(k_spacing))

    return kgrid.astype('int').tolist()
