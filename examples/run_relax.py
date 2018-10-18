import click
from click_spinner import spinner as cli_spinner
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.aiida_utils import load_dbenv_if_not_loaded, get_data_node
from auxiliary import example_param_set, create_structure_perturbed, set_kpoints, set_params_simple


@click.command()
@example_param_set
def main(potential_family, queue, code, computer):
    load_dbenv_if_not_loaded()
    from aiida.orm import WorkflowFactory, Code
    from aiida.work import submit, run

    # fetch the code to be used (tied to a computer)
    code = Code.get_from_string('{}@{}'.format(code, computer))

    # set the WorkChain you would like to call
    workflow = WorkflowFactory('vasp.relax')

    # organize options (needs a bit of special care)
    options = AttributeDict()
    options.account = ''
    options.qos = ''
    options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 20}
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
    inputs.code = code
    # set structure
    inputs.structure = create_structure_perturbed()
    # set parameters
    inputs.parameters = get_data_node('parameter', dict={
        'encut': 240, 'ismear': 0,
        'sigma': 0.1, 'system': 'test system',
    })  # yapf: disable
    # set k-point grid density
    inputs.kpoints = set_kpoints()
    # set potentials and their mapping
    inputs.potential_family = get_data_node('str', potential_family)
    inputs.potential_mapping = get_data_node('parameter', dict={'Si': 'Si'})
    # set options
    inputs.options = get_data_node('parameter', dict=options)
    # set settings
    inputs.settings = get_data_node('parameter', dict=settings)
    # set workchain related inputs
    # turn relaxation on
    inputs.perform = get_data_node('bool', True)
    inputs.convergence_on = get_data_node('bool', True)
    inputs.convergence_positions = get_data_node('float', 0.1)
    inputs.relax_parameters = get_data_node('parameter', dict={
        'ediffg': -0.0001
    })  # yapf: disable
    inputs.verbose = get_data_node('bool', True)
    # submit the requested workchain with the supplied inputs
    submit(workflow, **inputs)


if __name__ == '__main__':
    main()
