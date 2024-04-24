from time import sleep

import numpy as np
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import submit
from aiida.manage.configuration import load_profile
from aiida.orm import (
    Bool,
    Code,
    Float,
    Group,
    Int,
    QueryBuilder,
    Str,
    WorkChainNode,
    load_group,
)
from aiida.plugins import DataFactory, WorkflowFactory

load_profile()


def get_structure_SiC():
    """Set up SiC wurtzite cell

    wurtzite-type SiC
      1.0000000000
      3.0920000000   0.0000000000   0.0000000000
     -1.5460000000   2.6777505485   0.0000000000
      0.0000000000   0.0000000000   5.0730000000
    Si    C
        2     2
    Direct
      0.3333333333   0.6666666667   0.0000000000
      0.6666666667   0.3333333333   0.5000000000
      0.3333333333   0.6666666667   0.3758220000
      0.6666666667   0.3333333333   0.8758220000

    """

    StructureData = DataFactory('structure')
    a = 3.092
    c = 5.073
    lattice = [[a, 0, 0], [-a / 2, a / 2 * np.sqrt(3), 0], [0, 0, c]]
    structure = StructureData(cell=lattice)
    for pos_direct, symbol in zip(
        (
            [1.0 / 3, 2.0 / 3, 0],
            [2.0 / 3, 1.0 / 3, 0.5],
            [1.0 / 3, 2.0 / 3, 0.375822],
            [2.0 / 3, 1.0 / 3, 0.875822],
        ),
        ('Si', 'Si', 'C', 'C'),
    ):
        pos_cartesian = np.dot(pos_direct, lattice)
        structure.append_atom(position=pos_cartesian, symbols=symbol)
    return structure


def launch_aiida_relax_shape(structure, code_string, options, label):
    """Launch a relaxation of shape, but not volume."""

    Dict = DataFactory('dict')
    KpointsData = DataFactory('array.kpoints')

    # Set INCAR flags
    base_incar_dict = {
        'incar': {
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
    }

    # Set code, potentials to use and options
    base_config = {
        'code_string': code_string,
        'potential_family': 'PBE',
        'potential_mapping': {'Si': 'Si', 'C': 'C'},
        'options': options,
    }

    # Make sure some things are parsed
    base_parser_settings = {
        'add_energies': True,
        'add_forces': True,
        'add_stress': True,
    }

    # Set the provided code string
    code = Code.get_from_string(base_config['code_string'])

    # Which workchain to use
    workchain = WorkflowFactory('vasp.relax')

    # Construct the builder
    builder = workchain.get_builder()
    builder.code = code
    builder.parameters = Dict(dict=base_incar_dict)
    builder.structure = structure
    builder.settings = Dict(dict={'parser_settings': base_parser_settings})
    builder.potential_family = Str(base_config['potential_family'])
    builder.potential_mapping = Dict(dict=base_config['potential_mapping'])
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([6, 6, 4], offset=[0, 0, 0.5])
    builder.kpoints = kpoints
    builder.options = Dict(dict=base_config['options'])
    builder.metadata.label = label
    builder.metadata.description = label
    builder.clean_workdir = Bool(False)
    relax = AttributeDict()
    relax.perform = Bool(True)
    relax.force_cutoff = Float(1e-5)
    relax.positions = Bool(True)
    relax.shape = Bool(True)
    relax.volume = Bool(False)
    builder.relax = relax
    builder.verbose = Bool(True)

    # Submit calc
    node = submit(builder)
    return node


def launch_aiida_full_relax(structure, code_string, options, label):
    """Launch a relaxation of everything."""
    Dict = DataFactory('dict')
    KpointsData = DataFactory('array.kpoints')

    # Set INCAR flags
    base_incar_dict = {
        'incar': {
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
    }

    # Set code, potentials to use and options
    base_config = {
        'code_string': code_string,
        'potential_family': 'PBE.54',
        'potential_mapping': {'Si': 'Si', 'C': 'C'},
        'options': options,
    }

    # Make sure some things are parsed
    base_parser_settings = {
        'add_energies': True,
        'add_forces': True,
        'add_stress': True,
    }

    # Set the provided code string
    code = Code.get_from_string(base_config['code_string'])

    # Which workchain to use
    workchain = WorkflowFactory('vasp.relax')

    # Construct the builder
    builder = workchain.get_builder()
    builder.code = code
    builder.parameters = Dict(dict=base_incar_dict)
    builder.structure = structure
    builder.settings = Dict(dict={'parser_settings': base_parser_settings})
    builder.potential_family = Str(base_config['potential_family'])
    builder.potential_mapping = Dict(dict=base_config['potential_mapping'])
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([6, 6, 4], offset=[0, 0, 0.5])
    builder.kpoints = kpoints
    builder.options = Dict(dict=base_config['options'])
    builder.metadata.label = label
    builder.metadata.description = label
    builder.clean_workdir = Bool(False)
    relax = AttributeDict()
    relax.perform = Bool(True)
    relax.force_cutoff = Float(1e-5)
    relax.positions = Bool(True)
    relax.shape = Bool(True)
    relax.volume = Bool(True)
    relax.steps = Int(100)
    builder.relax = relax
    builder.verbose = Bool(True)

    # Submit the calc
    node = submit(builder)
    return node


def main(code_string, options, group_name, sleep_seconds=30):
    group = load_group(label=group_name)
    if not group.is_empty:
        # Make sure we do not already have members in the group
        print('The given group already have nodes. Please remove them or give a new group name.')
        exit(1)

    # Initial relaxation
    structure = get_structure_SiC()
    node_relax = launch_aiida_full_relax(structure, code_string, options, 'SiC VASP calc to relax volume')
    group.add_nodes(node_relax)
    print(f'Relaxing starting cell: {structure.cell}')

    while True:
        if node_relax.is_terminated:
            break
        print('Waiting for relaxation of starting cell to be done.')
        sleep(sleep_seconds)

    if node_relax.is_finished_ok:
        for strain, label in zip((0.99, 1.01), ('reduced', 'increased')):
            structure_clone = node_relax.outputs.relax.structure.clone()
            ase_structure = structure_clone.get_ase()
            # Make sure we also update atomic positions after scaling the cell
            ase_structure.set_cell(ase_structure.cell * strain ** (1.0 / 3.0), scale_atoms=True)
            StructureData = DataFactory('structure')
            structure = StructureData(ase=ase_structure)
            node = launch_aiida_relax_shape(
                structure,
                code_string,
                options,
                f'SiC VASP relax shape at {label} volume ({strain:f})',
            )
            group.add_nodes(node)
            print(f'Relaxation positions of scaled cell: {ase_structure.cell}')
    else:
        print('Relaxation failed.')
        exit(1)

    # Now make sure we wait until all calcs have been terminated
    while True:
        sleep(sleep_seconds)
        print('Waiting for all remaining relaxations to be done.')
        node_terminated = []
        for node in group.nodes:
            node_terminated.append(node.is_terminated)
        if all(node_terminated):
            break

    # Do one final check to make sure all calcs have a zero exit status
    for node in group.nodes:
        if not node.is_finished_ok:
            print(f'Node with id: {node.pk} has an exit status: {node.exit_status} and exit message: {node.exit_message}')
            exit(1)


def calc_bulk_modulus(group_name):
    stresses = []
    volumes = []
    for label in ('reduced', 'increased'):
        qb = QueryBuilder()
        qb.append(Group, filters={'label': group_name}, tag='group')
        qb.append(
            WorkChainNode,
            with_group='group',
            filters={'label': {'ilike': '%' + label + '%'}},
        )
        node = qb.first()[0]
        stresses.append(np.trace(node.outputs.stress.get_array('final')) / 3)
        volumes.append(np.linalg.det(node.inputs.structure.cell))

    d_s = stresses[1] - stresses[0]
    d_v = volumes[1] - volumes[0]
    v0 = (volumes[0] + volumes[1]) / 2.0
    bulk_modulus = -d_s / d_v * v0

    print(f'Bulk modulus: {bulk_modulus / 10.0} GPa')


if __name__ == '__main__':
    # Code_string is chosen from output of the list given by 'verdi code list'
    code_string = 'vasp@mycluster'

    # Set the options
    options = {
        'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 8},
        'account': '',
        'qos': '',
        'max_memory_kb': 2000000,
        'max_wallclock_seconds': 1800,
    }

    # Here we assume the group "Bulk_modulus_SiC_test" has already been created, but are empty
    group_name = 'Bulk_modulus_SiC_test'

    # Perform necessary calculations and put them in the created group
    main(code_string, options, group_name)

    # Calculate the bulk modulus by querying the group for the necessary input
    calc_bulk_modulus(group_name)
