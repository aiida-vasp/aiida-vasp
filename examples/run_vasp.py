import click
from click_spinner import spinner as cli_spinner
import numpy

from aiida_vasp.utils.aiida_utils import builder_interface


def get_data_cls(descriptor):
    load_dbenv_if_not_loaded()
    from aiida.orm import DataFactory
    return DataFactory(descriptor)


@click.group()
def run_example():
    """Run an example vasp calculation"""


def example_param_set(cmd_function):

    @click.option(
        '--pot-family', type=str, default='vasp-test', help='Name for a Potcar family to upload to (or to look for if --no-import).')
    @click.option(
        '--import-from',
        type=click.Path(),
        default='.',
        help='Folder to search for POTCARS to upload, . by default. Use --no-import to prevent uploading')
    @click.option('--no-import', is_flag=True, help='Do not try to import POTCARs, instead rely on the given family existing.')
    @click.option('--queue', type=str, default='', help='Name of the compute queue if your scheduler requires it')
    @click.argument('code', type=str)
    @click.argument('computer', type=str)
    def decorated_cmd_fn(*args, **kwargs):
        return cmd_function(*args, **kwargs)

    decorated_cmd_fn.__name__ = cmd_function.__name__
    decorated_cmd_fn.__doc__ = cmd_function.__doc__

    return decorated_cmd_fn


@run_example.command()
@example_param_set
def noncol(pot_family, import_from, queue, code, computer, no_import):
    load_dbenv_if_not_loaded()
    from aiida.orm import CalculationFactory, Code
    if not no_import:
        click.echo('importing POTCAR files...')
        with cli_spinner():
            import_pots(import_from, pot_family)
    pot_cls = get_data_cls('vasp.potcar')
    pots = {}
    pots['In'] = pot_cls.find_one(full_name='In_d', family=pot_family)
    pots['As'] = pot_cls.find_one(full_name='As', family=pot_family)

    vasp_calc = CalculationFactory('vasp.vasp')()
    vasp_calc.use_structure(create_structure_InAs())
    vasp_calc.use_kpoints(create_kpoints())
    vasp_calc.use_parameters(create_params_noncol())
    code = Code.get_from_string('{}@{}'.format(code, computer))
    vasp_calc.use_code(code)
    vasp_calc.use_potential(pots['In'], 'In')
    vasp_calc.use_potential(pots['As'], 'As')
    vasp_calc.set_computer(code.get_computer())
    vasp_calc.set_queue_name(queue)
    vasp_calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine': 20})
    vasp_calc.label = 'Test VASP run'
    vasp_calc.store_all()
    vasp_calc.submit()


@run_example.command()
@example_param_set
def simple(pot_family, import_from, queue, code, computer, no_import):
    load_dbenv_if_not_loaded()
    from aiida.orm import CalculationFactory, Code
    from aiida.work import submit
    if not no_import:
        click.echo('importing POTCAR files...')
        with cli_spinner():
            import_pots(import_from, pot_family)
    code = Code.get_from_string('{}@{}'.format(code, computer))

    calc_cls = CalculationFactory('vasp.vasp')
    if builder_interface(calc_cls):
        proc, inputs = simple_builder(pot_family, queue, code, calc_cls)
    else:
        proc, inputs = simple_inputs_template(pot_family, queue, code, calc_cls)

    submit(proc, **inputs)


def simple_inputs_template(pot_family, queue, code, calc_cls):
    """Submit simple example in AiiDA 0.10.0 - 0.11.4."""
    vasp_proc = calc_cls.process()
    inputs = vasp_proc.get_inputs_template()

    inputs.structure = create_structure_Si()
    inputs.kpoints = create_kpoints()
    inputs.parameters = create_params_simple()
    inputs.code = code
    inputs.potential = get_data_cls('vasp.potcar').get_potcars_from_structure(
        structure=builder.structure, family_name=pot_family, mapping={'Si': 'Si'})
    inputs._options.computer = builder.code.get_computer()
    inputs._options.queue_name = queue
    inputs._options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 8}
    inputs['_label'] = 'Example - simple VASP run'

    return vasp_proc, inputs


def simple_builder(pot_family, queue, code, calc_cls):
    """Submit simple example in AiiDA 0.12.0 and higher."""
    builder = calc_cls.get_builder()

    builder.structure = create_structure_Si()
    builder.kpoints = create_kpoints()
    builder.parameters = create_params_simple()
    builder.code = code
    builder.potential = get_data_cls('vasp.potcar').get_potcars_from_structure(
        structure=builder.structure, family_name=pot_family, mapping={'Si': 'Si'})
    builder.options.computer = builder.code.get_computer()
    builder.options.queue_name = queue
    builder.options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 8}
    builder.label = 'Example - simple VASP run'

    return calc_cls, builder


def load_dbenv_if_not_loaded():
    from aiida import load_dbenv, is_dbenv_loaded
    if not is_dbenv_loaded():
        load_dbenv()


def create_structure_InAs():
    structure_cls = get_data_cls('structure')
    structure = structure_cls(cell=numpy.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * 6.058)
    structure.append_atom(position=(0, 0, 0), symbols='In')
    structure.append_atom(position=(0.25, 0.25, 0.25), symbols='As')
    return structure


def create_structure_Si():
    structure_cls = get_data_cls('structure')
    alat = 5.4
    structure = structure_cls(cell=numpy.array([[.5, .5, 0], [.5, 0, .5], [0, .5, .5]]) * alat)
    structure.append_atom(position=numpy.array([.25, .25, .25]) * alat, symbols='Si')
    return structure


def create_kpoints():
    kpoints_cls = get_data_cls('array.kpoints')
    return kpoints_cls(kpoints_mesh=[8, 8, 8])


def create_params_noncol():
    param_cls = get_data_cls('parameter')
    return param_cls(
        dict={
            'SYSTEM': 'InAs',
            'EDIFF': 1e-5,
            'LORBIT': 11,
            'LSORBIT': '.True.',
            'GGA_COMPAT': '.False.',
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'GGA': 'PE',
            'ENCUT': '280.00 eV',
            'MAGMOM': '6*0.0',
            'NBANDS': 24,
        })


def create_params_simple():
    param_cls = get_data_cls('parameter')
    return param_cls(dict={'prec': 'NORMAL', 'encut': 200, 'ediff': 1e-8, 'ialgo': 38, 'ismear': 0, 'sigma': 0.1})


def import_pots(folder_path, family_name):
    pot_cls = get_data_cls('vasp.potcar')
    pot_cls.upload_potcar_family(folder_path, group_name=family_name, group_description='Test family', stop_if_existing=False)


if __name__ == '__main__':
    run_example()
