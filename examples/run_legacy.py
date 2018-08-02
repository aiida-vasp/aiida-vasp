"""This is an example of a legacy workflow."""

import click
from click_spinner import spinner as cli_spinner
import numpy

from aiida_vasp.utils.aiida_utils import builder_interface, load_dbenv_if_not_loaded
from auxilary import example_param_set, set_structure_InAs, set_structure_Si, set_kpoints, set_params_noncol, set_params_simple


@click.group()
def run_example():
    """Run an example vasp calculation"""


@run_example.command()
@example_param_set
def noncol(potcar_family, import_from, queue, code, computer, no_import):
    load_dbenv_if_not_loaded()
    from aiida.orm import CalculationFactory, Code
    from aiida.work import submit
    if not no_import:
        click.echo('importing POTCAR files...')
        with cli_spinner():
            import_pots(import_from, potcar_family)
    code = Code.get_from_string('{}@{}'.format(code, computer))
    calc_cls = CalculationFactory('vasp.vasp')

    if builder_interface(calc_cls):
        proc, inputs = noncol_builder(potcar_family, queue, code, calc_cls)
    else:
        proc, inputs = noncol_inputs_template(potcar_family, queue, code, calc_cls)

    submit(proc, **inputs)


def noncol_inputs_template(potcar_family, queue, code, calc_cls):
    """Submit noncol example in AiiDA 0.10.0 - 0.11.4."""
    vasp_proc = calc_cls.process()
    inputs = vasp_proc.get_inputs_template()

    inputs.structure = set_structure_InAs()
    inputs.kpoints = set_kpoints()
    inputs.parameters = set_params_noncol()
    inputs.code = code
    inputs.potential = get_data_class('vasp.potcar').get_potcars_from_structure(
        structure=inputs.structure, family_name=potcar_family, mapping={
            'In': 'In_d',
            'As': 'As'
        })
    inputs._options.computer = inputs.code.get_computer()
    inputs._options.queue_name = queue
    inputs._options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 8}
    inputs['_label'] = 'Example - noncol VASP run'

    return vasp_proc, inputs


def noncol_builder(potcar_family, queue, code, calc_cls):
    """Submit noncol example in AiiDA 0.12.0 and higher."""
    builder = calc_cls.get_builder()

    builder.structure = set_structure_InAs()
    builder.kpoints = set_kpoints()
    builder.parameters = set_params_noncol()
    builder.code = code
    builder.potential = get_data_class('vasp.potcar').get_potcars_from_structure(
        structure=builder.structure, family_name=potcar_family, mapping={
            'In': 'In_d',
            'As': 'As'
        })
    builder.options.computer = builder.code.get_computer()
    builder.options.queue_name = queue
    builder.options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 8}
    builder.label = 'Example - noncol VASP run'

    return calc_cls, builder


@run_example.command()
@example_param_set
def simple(potcar_family, import_from, queue, code, computer, no_import):
    load_dbenv_if_not_loaded()
    from aiida.orm import CalculationFactory, Code
    from aiida.work import submit
    if not no_import:
        click.echo('importing POTCAR files...')
        with cli_spinner():
            import_pots(import_from, potcar_family)
    code = Code.get_from_string('{}@{}'.format(code, computer))

    calc_cls = CalculationFactory('vasp.vasp')
    if builder_interface(calc_cls):
        proc, inputs = simple_builder(potcar_family, queue, code, calc_cls)
    else:
        proc, inputs = simple_inputs_template(potcar_family, queue, code, calc_cls)

    submit(proc, **inputs)


def simple_inputs_template(potcar_family, queue, code, calc_cls):
    """Submit simple example in AiiDA 0.10.0 - 0.11.4."""
    vasp_proc = calc_cls.process()
    inputs = vasp_proc.get_inputs_template()

    inputs.structure = set_structure_Si()
    inputs.kpoints = set_kpoints()
    inputs.parameters = set_params_simple()
    inputs.code = code
    inputs.potential = get_data_class('vasp.potcar').get_potcars_from_structure(
        structure=inputs.structure, family_name=potcar_family, mapping={'Si': 'Si'})
    inputs._options.computer = inputs.code.get_computer()
    inputs._options.queue_name = queue
    inputs._options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 8}
    inputs['_label'] = 'Example - simple VASP run'

    return vasp_proc, inputs


def simple_builder(potcar_family, queue, code, calc_cls):
    """Submit simple example in AiiDA 0.12.0 and higher."""
    builder = calc_cls.get_builder()

    builder.structure = set_structure_Si()
    builder.kpoints = set_kpoints()
    builder.parameters = set_params_simple()
    builder.code = code
    builder.potential = get_data_class('vasp.potcar').get_potcars_from_structure(
        structure=builder.structure, family_name=potcar_family, mapping={'Si': 'Si'})
    builder.options.computer = builder.code.get_computer()
    builder.options.queue_name = queue
    builder.options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 8}
    builder.label = 'Example - simple VASP run'

    return calc_cls, builder


def import_pots(folder_path, family_name):
    pot_cls = get_data_cls('vasp.potcar')
    pot_cls.upload_potcar_family(folder_path, group_name=family_name, group_description='Test family', stop_if_existing=False)


if __name__ == '__main__':
    run_example()
