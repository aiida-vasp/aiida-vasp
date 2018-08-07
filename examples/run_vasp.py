import click
import os
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.aiida_utils import load_dbenv_if_not_loaded, get_data_node
from auxiliary import example_param_set, set_structure_Si, set_kpoints, set_params_simple, set_params_simple_no_encut

os.system('verdi daemon restart')


@click.command()
@example_param_set
def main(potcar_family, queue, code, computer):
    load_dbenv_if_not_loaded()
    from aiida.orm import WorkflowFactory, Code
    from aiida.work import submit, run
    from aiida_vasp.utils.aiida_utils import get_data_class

    # set the code to be used, currently tied to a computer
    code = Code.get_from_string('{}@{}'.format(code, computer))

    # set the workchain you would like to call
    workchain = WorkflowFactory('vasp.verify')

    # set inputs
    inputs = AttributeDict()
    inputs.restart = AttributeDict()
    inputs.verify = AttributeDict()
    inputs.structure = set_structure_Si()
    inputs.kpoints = set_kpoints()
    #inputs.incar = set_params_simple_no_encut()
    inputs.incar = set_params_simple()
    inputs.code = code
    inputs.potcar_family = get_data_class('str')(potcar_family)
    inputs.potcar_mapping = get_data_node('parameter', dict={'Si': 'Si'})
    inputs.restart.max_iterations = get_data_class('int')(1)
    inputs.verify.max_iterations = get_data_class('int')(1)
    inputs.restart.clean_workdir = get_data_class('bool')(True)
    inputs.verify.clean_workdir = get_data_class('bool')(False)

    # set options
    options = AttributeDict()
    settings = AttributeDict()
    converge_settings = {'compress': False, 'displace': False}
    settings['converge'] = converge_settings
    if computer == 'unity':
        options.account = ''
        options.queue_name = ''
        options.qos = ''
        options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 20}
    elif computer == 'fram':
        options.account = 'nn9544k'
        options.queue_name = ''
        options.qos = 'devel'
        options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 16}
    inputs.options = get_data_node('parameter', dict=options)
    inputs.settings = get_data_node('parameter', dict=settings)

    # submit the requested workchain with the supplied inputs
    run(workchain, **inputs)


if __name__ == '__main__':
    main()
