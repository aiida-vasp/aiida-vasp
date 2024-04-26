"""
Call script to calculate the total energies for one volume of standard silicon.

This particular call script set up a standard calculation that execute a calculation for
the fcc silicon structure.
"""
from aiida import load_profile

# pylint: disable=too-many-arguments, invalid-name, import-outside-toplevel
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import submit
from aiida.orm import Bool, Code, Str
from aiida.plugins import DataFactory, WorkflowFactory

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
    workchain = WorkflowFactory('vasp.vasp')

    # And finally, we declare the options, settings and input containers
    settings = AttributeDict()
    inputs = AttributeDict()

    # Organize settings
    settings.parser_settings = {'add_dos': True}

    # Set inputs for the following WorkChain execution
    # Code
    inputs.code = Code.get_from_string(code_string)
    # Structure
    inputs.structure = structure
    # k-points grid density
    kpoints = kpoints_data()
    kpoints.set_kpoints_mesh(kmesh)
    inputs.kpoints = kpoints
    # Parameters
    inputs.parameters = dict_data(dict=incar)
    # Potential family and their mapping between element and potential type to use
    inputs.potential_family = Str(potential_family)
    inputs.potential_mapping = dict_data(dict=potential_mapping)
    # Options
    inputs.options = dict_data(dict=options)
    # Settings
    inputs.settings = dict_data(dict=settings)
    # Workchain related inputs, in this case, give more explicit output to report
    inputs.verbose = Bool(True)
    # Submit the workchain with the set inputs
    submit(workchain, **inputs)


if __name__ == '__main__':
    # Code_string is chosen among the list given by 'verdi code list'
    CODE_STRING = 'vasp@mycluster'

    # INCAR equivalent
    # Set input parameters
    INCAR = {'incar': {'encut': 240, 'ismear': -5, 'lorbit': 11}}

    # KPOINTS equivalent
    # Set kpoint mesh
    KMESH = [21, 21, 21]

    # POTCAR equivalent
    # Potential_family is chosen among the list given by
    # 'verdi data vasp-potcar listfamilies'
    POTENTIAL_FAMILY = 'PBE.54'
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
    OPTIONS.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 8}
    OPTIONS.queue_name = ''
    OPTIONS.max_wallclock_seconds = 3600
    OPTIONS.max_memory_kb = 2000000

    # POSCAR equivalent
    # Set the silicon structure
    STRUCTURE = get_structure('silicon_at_3_9')

    main(CODE_STRING, INCAR, KMESH, STRUCTURE, POTENTIAL_FAMILY, POTENTIAL_MAPPING, OPTIONS)
