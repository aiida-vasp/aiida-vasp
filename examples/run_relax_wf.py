import click
from click_spinner import spinner as cli_spinner
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.aiida_utils import load_dbenv_if_not_loaded, get_data_node
from run_vasp import example_param_set, create_structure_Si, create_kpoints, create_params_simple


@click.command()
@example_param_set
def main(pot_family, import_from, queue, code, computer, no_import):
    load_dbenv_if_not_loaded()
    from aiida.orm import WorkflowFactory, Code
    from aiida.work import submit
    from aiida.work.db_types import Str, Float, Bool

    if not no_import:
        click.echo('importing POTCAR files...')
        with cli_spinner():
            import_pots(import_from, pot_family)

    code = Code.get_from_string('{}@{}'.format(code, computer))
    workflow = WorkflowFactory('vasp.relax')

    inputs = AttributeDict()
    inputs.structure = create_structure_Si()
    inputs.kpoints = AttributeDict()
    inputs.kpoints.distance = Float(0.2)
    inputs.relax = AttributeDict()
    inputs.convergence = AttributeDict()
    inputs.convergence.shape = AttributeDict()
    inputs.convergence.on = Bool(True)
    inputs.convergence.positions = Float(0.1)
    inputs.restart = AttributeDict()
    inputs.code = code
    inputs.potcar_family = Str(pot_family)
    inputs.potcar_mapping = get_data_node('parameter', dict={'Si': 'Si'})
    options = AttributeDict()
    options.queue_name = queue
    options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 4}
    inputs.options = get_data_node('parameter', dict=options)

    submit(workflow, **inputs)


if __name__ == '__main__':
    main()
