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

    code = Code.get_from_string('{}@{}'.format(code, computer))
    workflow = WorkflowFactory('vasp.relax')

    inputs = AttributeDict()
    inputs.structure = create_structure_perturbed()
    inputs.kpoints = set_kpoints()
    inputs.relax = AttributeDict()
    inputs.relax.perform = get_data_node('bool', True)
    inputs.relax.convergence = AttributeDict()
    inputs.relax.convergence.shape = AttributeDict()
    inputs.relax.convergence.on = get_data_node('bool', True)
    inputs.relax.convergence.positions = get_data_node('float', 0.1)
    inputs.relax.incar = get_data_node('parameter', dict={
        'nsw': 1, 'ediffg': -0.0001, 'encut': 240, 'ismear': 0,
        'sigma': 0.1, 'system': 'test-case:test_relax_wf',
    })  # yapf: disable
    inputs.restart = AttributeDict()
    inputs.verify = AttributeDict()
    inputs.restart.max_iterations = get_data_node('int', 1)
    inputs.restart.clean_workdir = get_data_node('bool', False)
    inputs.verify.max_iterations = get_data_node('int', 1)
    inputs.verify.clean_workdir = get_data_node('bool', False)
    inputs.code = code
    inputs.potential_family = get_data_node('str', potential_family)
    inputs.potential_mapping = get_data_node('parameter', dict={'Si': 'Si'})
    options = AttributeDict()
    options.queue_name = queue
    options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 4}
    inputs.options = get_data_node('parameter', dict=options)

    run(workflow, **inputs)


if __name__ == '__main__':
    main()
