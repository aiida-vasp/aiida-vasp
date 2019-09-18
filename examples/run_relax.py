import numpy as np
from aiida.common.extendeddicts import AttributeDict
from aiida.orm import Code, Bool, Str, Float
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.engine import submit, run
from aiida import load_profile
load_profile()


def get_structure_Si():
    """Set up Si primitive cell

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

    StructureData = DataFactory('structure')
    alat = 5.431
    lattice = np.array([[.5, 0, .5], [.5, .5, 0], [0, .5, .5]]) * alat
    structure = StructureData(cell=lattice)
    for pos_direct in ([0.875, 0.875, 0.875],
                       [0.125, 0.125, 0.125]):
        pos_cartesian = np.dot(pos_direct, lattice)
        structure.append_atom(position=pos_cartesian, symbols='Si')
    return structure


def main(code_string, potential_family, resources):
    Dict = DataFactory('dict')

    # set the WorkChain you would like to call
    workflow = WorkflowFactory('vasp.relax')

    # organize options (needs a bit of special care)
    options = AttributeDict()
    options.account = ''
    options.qos = ''
    options.resources = resources
    options.queue_name = ''
    options.max_wallclock_seconds = 3600

    # organize settings
    settings = AttributeDict()
    # the workchains should configure the required parser settings on the fly
    parser_settings = {}
    settings.parser_settings = parser_settings

    # set inputs for the following WorkChain execution

    inputs = AttributeDict()
    # set code
    inputs.code = Code.get_from_string(code_string)
    # set structure
    inputs.structure = get_structure_Si()
    # set parameters
    inputs.parameters = Dict(dict={'encut': 240,
                                   'ismear': 0,
                                   'sigma': 0.1,
                                   'system': 'test system'})
    # set k-point grid density
    KpointsData = DataFactory("array.kpoints")
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([9, 9, 9])
    inputs.kpoints = kpoints
    # set potentials and their mapping
    inputs.potential_family = Str(potential_family)
    inputs.potential_mapping = Dict(dict={'Si': 'Si'})
    # set options
    inputs.options = Dict(dict=options)
    # set settings
    inputs.settings = Dict(dict=settings)
    # set workchain related inputs
    # turn relaxation on
    inputs.relax = Bool(True)
    inputs.force_cutoff = Float(0.01)
    inputs.convergence_on = Bool(False)
    inputs.convergence_positions = Float(0.1)
    inputs.relax_parameters = Dict(dict={'ediffg': -0.0001,
                                         'ibrion': 2,
                                         'nsw': 10,
                                         'isif': 3})
    inputs.verbose = Bool(True)
    # submit the requested workchain with the supplied inputs
    submit(workflow, **inputs)


if __name__ == '__main__':
    # code_string is chosen among the list given by 'verdi code list'
    code_string = 'vasp544mpi@mycomputer'

    # potential_family is chosen among the list given by
    # 'verdi data vasp-potcar listfamilies'
    potential_family = 'PBE.54'

    # metadata.options.resources
    # See https://aiida.readthedocs.io/projects/aiida-core/en/latest/scheduler/index.html
    resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 20}
    # resources = {'parallel_env': 'mpi*', 'tot_num_mpiprocs': 12}

    main(code_string, potential_family, resources)
