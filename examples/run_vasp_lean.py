import click
import os
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.aiida_utils import load_dbenv_if_not_loaded, get_data_node
from auxiliary import example_param_set, set_structure_si, set_kpoints, set_params_simple, set_params_simple_no_encut

os.system('verdi daemon restart')


@click.command()
@example_param_set
def main(potential_family, queue, code, computer):
    load_dbenv_if_not_loaded()
    from aiida.orm import WorkflowFactory, Code
    from aiida.work import submit, run
    from aiida_vasp.utils.aiida_utils import get_data_class

    # set the code to be used, currently tied to a computer
    code = Code.get_from_string('{}@{}'.format(code, computer))

    # set the workchain you would like to call
    workchain = WorkflowFactory('vasp.verify')

    # organize options (needs a bit of special care)
    options = AttributeDict()
    options.account = ''
    options.qos = ''
    options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 20}
    options.queue_name = ''
    options.max_wallclock_seconds = 3600

    # organize selfettings
    settings = AttributeDict()
    parser_settings = {'output_params': ['total_energies', 'maximum_force']}
    settings.parser_settings = parser_settings

    # set inputs for the following WorkChain execution

    inputs = AttributeDict()
    # set code
    inputs.code = code
    # set structure
    inputs.structure = set_structure_si()
    # set k-points grid density
    inputs.kpoints = set_kpoints()
    # set parameters
    inputs.parameters = set_params_simple()
    # set potentials and their mapping
    inputs.potential_family = get_data_class('str')(potential_family)
    inputs.potential_mapping = get_data_node('parameter', dict={'Si': 'Si'})
    # set options
    inputs.options = get_data_node('parameter', dict=options)
    # set settings
    inputs.settings = get_data_node('parameter', dict=settings)
    # set workchain related inputs
    inputs.verbose = get_data_node('bool', True)
    # submit the requested workchain with the supplied inputs
    run(workchain, **inputs)


if __name__ == '__main__':
    main()
