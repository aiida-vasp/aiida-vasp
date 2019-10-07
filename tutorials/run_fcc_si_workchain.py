"""
Call script to calculate the total energies for different volumes of the silicon structure.

This particular call script calls the workchain that handles the execution of the
calculation of each structure. In this script we only supply a dictionary of different
structures as an input to the workchain, in addition to the standard inputs required by
the AiiDA-VASP workchain stack.
"""
# pylint: disable=too-many-arguments
import numpy as np
from aiida.common.extendeddicts import AttributeDict
from aiida.orm import Code, Bool, Str
from aiida.plugins import DataFactory
from aiida.engine import submit
from aiida import load_profile
from eos import EosWorkChain
load_profile()


def get_structure(lattice_constant):
    """
    Set up Si primitive cell

    fcc Si:
       3.9
       0.5000000000000000    0.5000000000000000    0.0000000000000000
       0.0000000000000000    0.5000000000000000    0.5000000000000000
       0.5000000000000000    0.0000000000000000    0.5000000000000000
    Si
       1
    Cartesian
    0.0000000000000000  0.0000000000000000  0.0000000000000000

    """

    structure_data = DataFactory('structure')
    alat = lattice_constant
    lattice = np.array([[.5, .5, 0], [0, .5, .5], [.5, 0, .5]]) * alat
    structure = structure_data(cell=lattice)
    positions = [[0.0, 0.0, 0.0]]
    for pos_direct in positions:
        pos_cartesian = np.dot(pos_direct, lattice)
        structure.append_atom(position=pos_cartesian, symbols='Si')
    return structure


def get_structures(lattice_constants):
    """Build a dictionary of structures using different lattice constants."""
    structures = {}
    for lattice_constant in lattice_constants:
        # use the lattice constant as a key
        structures[str(lattice_constant).replace('.', '_')] = get_structure(lattice_constant)
    return structures


def main(code_string, incar, kmesh, structures, potential_family, potential_mapping, options):
    """Main method to setup the calculation."""

    # First, we need to fetch the AiiDA datatypes which will
    # house the inputs to our calculation
    dict_data = DataFactory('dict')
    kpoints_data = DataFactory('array.kpoints')

    # Then, we set the workchain you would like to call
    #workchain = WorkflowFactory('vasp.verify')
    workchain = EosWorkChain

    # And finally, we declare the options, settings and input containers
    settings = AttributeDict()
    inputs = AttributeDict()

    # organize settings
    settings.parser_settings = {'output_params': ['total_energies', 'maximum_force']}

    # set inputs for the following WorkChain execution
    # set code
    inputs.code = Code.get_from_string(code_string)
    # set structures
    inputs.structures = structures
    # set k-points grid density
    kpoints = kpoints_data()
    kpoints.set_kpoints_mesh(kmesh)
    inputs.kpoints = kpoints
    # set parameters
    inputs.parameters = dict_data(dict=incar)
    # set potentials and their mapping
    inputs.potential_family = Str(potential_family)
    inputs.potential_mapping = dict_data(dict=potential_mapping)
    # set options
    inputs.options = dict_data(dict=options)
    # set settings
    inputs.settings = dict_data(dict=settings)
    # set workchain related inputs, in this case, give more explicit output to report
    inputs.verbose = Bool(True)
    # submit the requested workchain with the supplied inputs
    submit(workchain, **inputs)


if __name__ == '__main__':
    # Code_string is chosen among the list given by 'verdi code list'
    CODE_STRING = 'vasp@mycluster'

    # INCAR equivalent
    # Set input parameters
    INCAR = {'istart': 0, 'icharg': 2, 'encut': 240, 'ismear': 0, 'sigma': 0.1}

    # KPOINTS equivalent
    # Set kpoint mesh
    KMESH = [11, 11, 11]

    # POTCAR equivalent
    # Potential_family is chosen among the list given by
    # 'verdi data vasp-potcar listfamilies'
    POTENTIAL_FAMILY = 'pbe'
    # The potential mapping selects which potential to use, here we use the standard
    # for silicon, this could for instance be {'Si': 'Si_GW'} to use the GW ready
    # potential instead
    POTENTIAL_MAPPING = {'Si': 'Si'}

    # jobfile equivalent
    # In options, we typically set scheduler options.
    # See https://aiida.readthedocs.io/projects/aiida-core/en/latest/scheduler/index.html
    # AttributeDict is just a special dictionary with the extra benefit that
    # you can set and get the key contents with mydict.mykey, instead of mydict['mykey']
    OPTIONS = AttributeDict()
    OPTIONS.account = 'nn9995k'
    OPTIONS.qos = ''
    OPTIONS.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 16}
    OPTIONS.queue_name = ''
    OPTIONS.max_wallclock_seconds = 3600
    OPTIONS.max_memory_kb = 1024000

    # POSCAR equivalent
    # Set the silicon structure
    LATTICE_CONSTANTS = [3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3]
    STRUCTURES = get_structures(LATTICE_CONSTANTS)

    main(CODE_STRING, INCAR, KMESH, STRUCTURES, POTENTIAL_FAMILY, POTENTIAL_MAPPING, OPTIONS)
