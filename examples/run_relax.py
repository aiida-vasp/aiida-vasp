"""
An example call script that performs a relaxation of a structure.

Performs a relaxation of the standard silicon structure.
"""
# pylint: disable=too-many-arguments
import numpy as np
from aiida.common.extendeddicts import AttributeDict
from aiida.orm import Code
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.engine import submit
from aiida import load_profile
load_profile()


def get_structure():
    """
    Set up Si primitive cell

    Si
       5.431
         0.0000000000000000    0.5000000000000000    0.5000000000000000
         0.4900000000000000    0.0000000000000000    0.4800000000000000
         0.5000000000000000    0.5000000000000000    0.0100000000000000
    Si
       2
    Direct
      0.8750000000000000  0.8750000000000000  0.8750000000000000
      0.1250000000000000  0.1250000000000000  0.1250000000000000

    """

    structure_data = DataFactory('structure')
    alat = 5.431
    lattice = np.array([[.5, 0, .5], [.5, .5, 0], [0, .5, .5]]) * alat
    structure = structure_data(cell=lattice)
    for pos_direct in ([0.875, 0.875, 0.875], [0.125, 0.125, 0.125]):
        pos_cartesian = np.dot(pos_direct, lattice)
        structure.append_atom(position=pos_cartesian, symbols='Si')
    return structure


def main(code_string, incar, kmesh, structure, potential_family, potential_mapping, options):
    """Main method to setup the calculation."""

    # We set the workchain you would like to call
    workchain = WorkflowFactory('vasp.relax')

    # And finally, we declare the options, settings and input containers
    settings = AttributeDict()
    inputs = AttributeDict()

    # Organize settings
    settings.parser_settings = {}

    # Set inputs for the following WorkChain execution
    # Set code
    inputs.code = Code.get_from_string(code_string)
    # Set structure
    inputs.structure = structure
    # Set k-points grid density
    kpoints = DataFactory('array.kpoints')()
    kpoints.set_kpoints_mesh(kmesh)
    inputs.kpoints = kpoints
    # Set parameters
    inputs.parameters = DataFactory('dict')(dict=incar)
    # Set potentials and their mapping
    inputs.potential_family = DataFactory('str')(potential_family)
    inputs.potential_mapping = DataFactory('dict')(dict=potential_mapping)
    # Set options
    inputs.options = DataFactory('dict')(dict=options)
    # Set settings
    inputs.settings = DataFactory('dict')(dict=settings)
    # Set workchain related inputs, in this case, give more explicit output to report
    inputs.verbose = DataFactory('bool')(True)

    # Relaxation related parameters that is passed to the relax workchain
    relax = AttributeDict()
    # Turn on relaxation
    relax.perform = DataFactory('bool')(True)
    # Select relaxation algorithm
    relax.algo = DataFactory('str')('cg')
    # Set force cutoff limit (EDIFFG, but no sign needed)
    relax.force_cutoff = DataFactory('float')(0.01)
    # Turn on relaxation of positions (strictly not needed as the default is on)
    # The three next parameters correspond to the well known ISIF=3 setting
    relax.positions = DataFactory('bool')(True)
    # Turn on relaxation of the cell shape (defaults to False)
    relax.shape = DataFactory('bool')(True)
    # Turn on relaxation of the volume (defaults to False)
    relax.volume = DataFactory('bool')(True)
    # Set maximum number of ionic steps
    relax.steps = DataFactory('int')(100)
    # Set the relaxation parameters on the inputs
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
    # Set input parameters
    INCAR = {'vasp': {'encut': 240, 'ismear': 0, 'sigma': 0.1, 'system': 'test system'}}

    # KPOINTS equivalent
    # Set kpoint mesh
    KMESH = [9, 9, 9]

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
    OPTIONS.max_wallclock_seconds = 3600
    OPTIONS.max_memory_kb = 1024000

    main(CODE_STRING, INCAR, KMESH, STRUCTURE, POTENTIAL_FAMILY, POTENTIAL_MAPPING, OPTIONS)
