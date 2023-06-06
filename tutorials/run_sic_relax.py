import numpy as np

from aiida.common.extendeddicts import AttributeDict
from aiida.engine import submit
from aiida.manage.configuration import load_profile
from aiida.orm import Bool, Code, Float, Int, Str
from aiida.plugins import DataFactory, WorkflowFactory

load_profile()


def launch_aiida(structure, code_string, options, potential_family, potential_mapping, label='SiC VASP calculation'):
    Dict = DataFactory('dict')
    KpointsData = DataFactory('array.kpoints')

    incar_dict = {
        'PREC': 'Accurate',
        'EDIFF': 1e-8,
        'NELMIN': 5,
        'NELM': 100,
        'ENCUT': 500,
        'IALGO': 38,
        'ISMEAR': 0,
        'SIGMA': 0.01,
        'LREAL': False,
        'LCHARG': False,
        'LWAVE': False,
    }

    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([6, 6, 4], offset=[0, 0, 0.5])

    parser_settings = {'add_energies': True, 'add_forces': True, 'add_stress': True}

    code = Code.get_from_string(code_string)
    Workflow = WorkflowFactory('vasp.relax')
    builder = Workflow.get_builder()
    builder.code = code
    builder.parameters = Dict(dict={'incar': incar_dict})
    builder.structure = structure
    builder.settings = Dict(dict={'parser_settings': parser_settings})
    builder.potential_family = Str(potential_family)
    builder.potential_mapping = Dict(dict=potential_mapping)
    builder.kpoints = kpoints
    builder.options = Dict(dict=options)
    builder.metadata.label = label
    builder.metadata.description = label
    relax = AttributeDict()
    relax.perform = Bool(True)  # Turn on relaxation of the structure
    relax.force_cutoff = Float(1e-5)  # Relax force cutoff
    relax.steps = Int(10)  # Relax number of ionic steps cutoff
    relax.positions = Bool(True)  # Relax atomic positions
    relax.shape = Bool(True)  # Relax cell shape (alpha, beta, gamma)
    relax.volume = Bool(True)  # Relax volume
    builder.relax = relax
    builder.verbose = Bool(True)
    builder.clean_workdir = Bool(False)
    node = submit(builder)
    return node


def get_structure_SiC():
    """Set up SiC cell

    Si C
       1.0
         3.0920072935808083    0.0000000000000000    0.0000000000000000
        -1.5460036467904041    2.6777568649277486    0.0000000000000000
         0.0000000000000000    0.0000000000000000    5.0733470000000001
     Si C
       2   2
    Direct
       0.3333333333333333  0.6666666666666665  0.4995889999999998
       0.6666666666666667  0.3333333333333333  0.9995889999999998
       0.3333333333333333  0.6666666666666665  0.8754109999999998
       0.6666666666666667  0.3333333333333333  0.3754109999999997

    """

    StructureData = DataFactory('structure')
    a = 3.092
    c = 5.073
    lattice = [[a, 0, 0], [-a / 2, a / 2 * np.sqrt(3), 0], [0, 0, c]]
    structure = StructureData(cell=lattice)
    for pos_direct, symbol in zip(([1. / 3, 2. / 3, 0], [2. / 3, 1. / 3, 0.5], [1. / 3, 2. / 3, 0.375822], [2. / 3, 1. / 3, 0.875822]),
                                  ('Si', 'Si', 'C', 'C')):
        pos_cartesian = np.dot(pos_direct, lattice)
        structure.append_atom(position=pos_cartesian, symbols=symbol)
    return structure


def main(code_string, options, potential_family, potential_mapping):
    structure = get_structure_SiC()
    launch_aiida(structure, code_string, options, potential_family, potential_mapping)


if __name__ == '__main__':
    code_string = 'vasp@mycluster'
    options = {
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 8
        },
        'account': '',
        'max_memory_kb': 2000000,
        'max_wallclock_seconds': 3600
    }
    potential_family = 'PBE.54'
    potential_mapping = {'Si': 'Si', 'C': 'C'}
    main(code_string, options, potential_family, potential_mapping)
