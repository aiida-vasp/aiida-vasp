"""Contains auxiliary function in order to be able to run the example workchain and workflows."""
import click
import numpy

from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node


def example_param_set(cmd_function):
    """Define the command line options."""

    @click.option(
        '--potcar-family', type=str, default='vasp-test', help='Name for a Potcar family to upload to (or to look for if --no-import).')
    @click.option('--queue', type=str, default='', help='Name of the compute queue if your scheduler requires it')
    @click.argument('code', type=str)
    @click.argument('computer', type=str)
    def decorated_cmd_fn(*args, **kwargs):
        return cmd_function(*args, **kwargs)

    decorated_cmd_fn.__name__ = cmd_function.__name__
    decorated_cmd_fn.__doc__ = cmd_function.__doc__

    return decorated_cmd_fn


def create_structure_perturbed():
    """
    Create a perturbed structure (example taken from the VASP wiki).


    `Link to the wiki page <http://cms.mpi.univie.ac.at/wiki/index.php/Cd_Si_relaxation>`_
    """

    alat = 5.5
    structure = get_data_node('structure', cell=numpy.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * alat)
    structure.append_atom(position=numpy.array([-0.125, -0.125, -0.125]) * alat, symbols='Si')
    structure.append_atom(position=numpy.array([0.125, 0.125, 0.130]) * alat, symbols='Si')
    return structure


def set_structure_InAs():
    """Set up a simple InAs structure."""

    structure_cls = get_data_class('structure')
    structure = structure_cls(cell=numpy.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * 6.058)
    structure.append_atom(position=(0, 0, 0), symbols='In')
    structure.append_atom(position=(0.25, 0.25, 0.25), symbols='As')
    return structure


def set_structure_Si():
    """Set up a simple Si structure."""

    structure_cls = get_data_class('structure')
    alat = 5.4
    structure = structure_cls(cell=numpy.array([[.5, 0, .5], [.5, .5, 0], [0, .5, .5]]) * alat)
    structure.append_atom(position=numpy.array([.25, .25, .25]) * alat, symbols='Si')
    return structure


def set_kpoints():
    """Set a simple kpoint sampling."""

    kpoints_cls = get_data_class('array.kpoints')
    return kpoints_cls(kpoints_mesh=[8, 8, 8])


def set_params_noncol():
    """Set INCAR parameters."""

    param_cls = get_data_class('parameter')
    return param_cls(
        dict={
            'SYSTEM': 'InAs',
            'EDIFF': 1e-5,
            'LORBIT': 11,
            'LSORBIT': '.True.',
            'GGA_COMPAT': '.False.',
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'GGA': 'PE',
            'ENCUT': '280.00 eV',
            'MAGMOM': '6*0.0',
            'NBANDS': 24,
        })


def set_params_simple():
    """Set INCAR parameters."""

    param_cls = get_data_class('parameter')
    return param_cls(dict={'prec': 'NORMAL', 'encut': 200, 'ediff': 1E-4, 'ialgo': 38, 'ismear': -5, 'sigma': 0.1})


def set_params_simple_no_encut():
    """Set INCAR parameters."""

    param_cls = get_data_class('parameter')
    return param_cls(dict={'prec': 'NORMAL', 'ediff': 1E-4, 'ialgo': 38, 'ismear': 0, 'sigma': 0.1})
