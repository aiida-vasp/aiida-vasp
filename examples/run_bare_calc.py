"""This is an example of calculation like that in tutorial for QE."""

import click
from aiida_vasp.utils.aiida_utils import load_dbenv_if_not_loaded
from auxiliary import (example_param_set, set_structure_si,
                       set_kpoints, set_params_simple)


@click.command()
@example_param_set
def main(potential_family, queue, code, computer):
    load_dbenv_if_not_loaded()
    from aiida.orm import Code, DataFactory
    from aiida.work import submit

    code = Code.get_from_string('{}@{}'.format(code, computer))
    builder = code.get_builder()
    builder.structure = set_structure_si()
    builder.kpoints = set_kpoints()
    builder.parameters = set_params_simple()
    builder.potential = DataFactory('vasp.potcar').get_potcars_from_structure(
        structure=builder.structure,
        family_name=potential_family,
        mapping={'Si': 'Si'})
    builder.options.queue_name = queue
    builder.options.resources = {'num_machines': 1,
                                 'num_mpiprocs_per_machine': 20}
    builder.label = "simple VASP run"
    builder.description = "Example - simple VASP run"

    settings_dict = {'parser_settings': {'add_forces': True}}
    builder.settings = DataFactory('parameter')(dict=settings_dict)

    submit(builder)


if __name__ == '__main__':
    main()
