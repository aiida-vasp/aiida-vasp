import numpy
import click
from click_spinner import spinner as cli_spinner
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.aiida_utils import load_dbenv_if_not_loaded, get_data_node
from run_vasp import example_param_set, create_structure_Si, create_kpoints, create_params_simple


def create_structure_perturbed():
    """
    Create a perturbed structure (example taken from the VASP wiki).


    `Link to the wiki page <http://cms.mpi.univie.ac.at/wiki/index.php/Cd_Si_relaxation>`_
    """

    alat = 5.5
    structure = get_data_node('structure', cell=numpy.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * alat)
    structure.append_atom(position=numpy.array([-0.125, -0.125, -0.125]) * alat, symbols='Si')
    structure.append_atom(position=numpy.array([0.125, 0.125, 0.130]) * alat, symbols='Si')
    return structure


@click.command()
@example_param_set
def main(pot_family, import_from, queue, code, computer, no_import):
    load_dbenv_if_not_loaded()
    from aiida.orm import WorkflowFactory, Code
    from aiida.work import submit

    if not no_import:
        click.echo('importing POTCAR files...')
        with cli_spinner():
            import_pots(import_from, pot_family)

    code = Code.get_from_string('{}@{}'.format(code, computer))
    workflow = WorkflowFactory('vasp.relax')

    inputs = AttributeDict()
    inputs.structure = create_structure_perturbed()
    inputs.kpoints = AttributeDict()
    inputs.kpoints.distance = get_data_node('float', 0.2)
    inputs.relax = AttributeDict()
    inputs.convergence = AttributeDict()
    inputs.convergence.shape = AttributeDict()
    inputs.convergence.on = get_data_node('bool', True)
    inputs.convergence.positions = get_data_node('float', 0.1)
    inputs.incar_add = get_data_node('parameter', dict={
        'nsw': 1, 'ediffg': -0.0001, 'encut': 240, 'ismear': 0, 'sigma': 0.1, 'system': 'test-case:test_relax_wf',
    })  # yapf: disable
    inputs.restart = AttributeDict()
    inputs.code = code
    inputs.potcar_family = get_data_node('str', pot_family)
    inputs.potcar_mapping = get_data_node('parameter', dict={'Si': 'Si'})
    options = AttributeDict()
    options.queue_name = queue
    options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 4}
    inputs.options = get_data_node('parameter', dict=options)

    submit(workflow, **inputs)


if __name__ == '__main__':
    main()
