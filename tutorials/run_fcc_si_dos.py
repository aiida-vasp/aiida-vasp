"""
Call script to calculate the total energies for one volume of standard silicon.

This particular call script set up a standard calculation that execute a calculation for
the fcc silicon structure.
"""
# pylint: disable=too-many-arguments, invalid-name, import-outside-toplevel
from aiida.common.extendeddicts import AttributeDict
from aiida.orm import Code, Bool, Str
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.engine import submit
from aiida import load_profile
load_profile()


def get_structure(label):
    from aiida.orm import QueryBuilder
    qb = QueryBuilder()
    qb.append(DataFactory('structure'), filters={'label': {'==': label}}, tag='structure')
    # Pick any structure with this label, here, just the first
    return qb.all()[0][0]


def main(code_string, incar, kmesh, structure, potential_family, potential_mapping, options):
    """Main method to setup the calculation."""

    # First, we need to fetch the AiiDA datatypes which will
    # house the inputs to our calculation
    dict_data = DataFactory('dict')
    kpoints_data = DataFactory('array.kpoints')

    # Then, we set the workchain you would like to call
    workchain = WorkflowFactory('vasp.verify')

    # And finally, we declare the options, settings and input containers
    settings = AttributeDict()
    inputs = AttributeDict()

    # organize settings
    settings.parser_settings = {'add_dos': True}

    # set inputs for the following WorkChain execution
    # set code
    inputs.code = Code.get_from_string(code_string)
    # set structure
    inputs.structure = structure
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
    INCAR = {'encut': 240, 'ismear': -5, 'lorbit': 11}

    # KPOINTS equivalent
    # Set kpoint mesh
    KMESH = [21, 21, 21]

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
    OPTIONS.account = ''
    OPTIONS.qos = ''
    OPTIONS.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 16}
    OPTIONS.queue_name = ''
    OPTIONS.max_wallclock_seconds = 3600
    OPTIONS.max_memory_kb = 1024000

    # POSCAR equivalent
    # Set the silicon structure
    STRUCTURE = get_structure('silicon_at_3_9')

    main(CODE_STRING, INCAR, KMESH, STRUCTURE, POTENTIAL_FAMILY, POTENTIAL_MAPPING, OPTIONS)
