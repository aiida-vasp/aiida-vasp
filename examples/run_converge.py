"""
An example call script that performs a single static VASP calculation.

Performs a self consistent electron convergence run using the standard silicon structure.
"""
# pylint: disable=too-many-arguments
import numpy as np
from aiida import load_profile
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import submit
from aiida.orm import Bool, Code, Str
from aiida.plugins import DataFactory, WorkflowFactory

load_profile()


def get_structure():
    """
    Set up Si primitive cell

    Si
       5.431
         0.0000000000000000    0.5000000000000000    0.5000000000000000
         0.5000000000000000    0.0000000000000000    0.5000000000000000
         0.5000000000000000    0.5000000000000000    0.0000000000000000
    Si
       2
    Direct
      0.8750000000000000  0.8750000000000000  0.8750000000000000
      0.1250000000000000  0.1250000000000000  0.1250000000000000

    """

    structure_data = DataFactory('core.structure')
    alat = 5.431
    lattice = np.array([[0.5, 0, 0.5], [0.5, 0.5, 0], [0, 0.5, 0.5]]) * alat
    structure = structure_data(cell=lattice)
    for pos_direct in ([0.875, 0.875, 0.875], [0.125, 0.125, 0.125]):
        pos_cartesian = np.dot(pos_direct, lattice)
        structure.append_atom(position=pos_cartesian, symbols='Si')
    return structure


def main(code_string, incar, kmesh, structure, potential_family, potential_mapping, options):
    """Main method to setup the calculation."""

    # First, we need to fetch the AiiDA datatypes which will
    # house the inputs to our calculation
    dict_data = DataFactory('core.dict')
    kpoints_data = DataFactory('core.array.kpoints')

    # Then, we set the workchain you would like to call
    workchain = WorkflowFactory('vasp.converge')

    # And finally, we declare the options, settings and input containers
    settings = AttributeDict()
    inputs = AttributeDict()

    # Organize settings
    settings.parser_settings = {'output_params': ['total_energies', 'maximum_force']}

    # Set inputs for the following WorkChain execution
    # Set code
    inputs.code = Code.get_from_string(code_string)
    # Set structure
    inputs.structure = structure
    # Set k-points grid density
    if kmesh:
        # Only set it if kmesh is supplied, otherwise run convergence
        # tests
        kpoints = kpoints_data()
        kpoints.set_kpoints_mesh(kmesh)
        inputs.kpoints = kpoints
    # Set parameters
    inputs.parameters = dict_data(dict=incar)
    # Set potentials and their mapping
    inputs.potential_family = Str(potential_family)
    inputs.potential_mapping = dict_data(dict=potential_mapping)
    # Set options
    inputs.options = dict_data(dict=options)
    # Set settings
    inputs.settings = dict_data(dict=settings)
    # Set workchain related inputs, in this case, give more explicit output to repor
    inputs.verbose = Bool(True)

    # Convergence and relaxation related parameters that is passed to the convergence
    # and relaxation workchain, respectively
    # Turn of final relaxation (after convergence tests)
    relax = AttributeDict()
    relax.perform = Bool(False)
    inputs.relax = relax
    # Submit the requested workchain with the supplied inputs
    submit(workchain, **inputs)


if __name__ == '__main__':
    # Code_string is chosen among the list given by 'verdi code list'
    CODE_STRING = 'vasp@mycluster'

    # POSCAR equivalent
    # Set the silicon structure
    STRUCTURE = get_structure()

    # INCAR equivalent
    # Set input parameters (make sure we do not set ENCUT in order to run convergence tests
    # on the plane wave cutoff)
    INCAR = {'incar': {'prec': 'NORMAL', 'ediff': 1e-4, 'ialgo': 38, 'ismear': -5, 'sigma': 0.1}}

    # KPOINTS equivalent
    # Set kpoint mesh
    KMESH = []

    # POTCAR equivalent
    # Potential_family is chosen among the list given by
    # 'verdi data vasp-potcar listfamilies'
    POTENTIAL_FAMILY = 'pbe'
    # The potential mapping selects which potential to use, here we use the standard
    # for silicon, this could for instance be {'Si': 'Si_GW'} to use the GW ready
    # potential instead
    POTENTIAL_MAPPING = {'Si': 'Si'}

    # Jobfile equivalent
    # In options, we typically set scheduler options.
    # See https://aiida.readthedocs.io/projects/aiida-core/en/latest/scheduler/index.html
    # AttributeDict is just a special dictionary with the extra benefit that
    # you can set and get the key contents with mydict.mykey, instead of mydict['mykey']
    OPTIONS = AttributeDict()
    OPTIONS.account = ''
    OPTIONS.qos = ''
    OPTIONS.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 1}
    OPTIONS.queue_name = ''
    OPTIONS.max_wallclock_seconds = 86400
    OPTIONS.max_memory_kb = 9000000

    main(CODE_STRING, INCAR, KMESH, STRUCTURE, POTENTIAL_FAMILY, POTENTIAL_MAPPING, OPTIONS)
