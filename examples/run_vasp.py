import click
import numpy


@click.group()
def run_example():
    """Run an example vasp calculation"""


def example_param_set(cmd_function):

    @click.option('--paw-family', type=str, default='vasp-test')
    @click.option('--import-from', type=click.Path(), default='.')
    @click.option('--queue', type=str, default='')
    @click.argument('code', type=str)
    @click.argument('computer', type=str)
    def decorated_cmd_fn(*args, **kwargs):
        return cmd_function(*args, **kwargs)

    decorated_cmd_fn.__name__ = cmd_function.__name__
    decorated_cmd_fn.__doc__ = cmd_function.__doc__

    return decorated_cmd_fn


@run_example.command()
@example_param_set
def noncol(paw_family, import_from, queue, code, computer):
    load_dbenv_if_not_loaded()
    from aiida.orm import CalculationFactory, Code
    if import_from:
        import_paws(import_from, paw_family)
    paw_in, paw_as = get_paws(paw_family, 'In', 'As')

    vasp_calc = CalculationFactory('vasp.vasp')()
    vasp_calc.use_structure(create_structure_InAs())
    vasp_calc.use_kpoints(create_kpoints())
    vasp_calc.use_parameters(create_params())
    code = Code.get_from_string('{}@{}'.format(code, computer))
    vasp_calc.use_code(code)
    vasp_calc.use_paw(paw_in, 'In')
    vasp_calc.use_paw(paw_as, 'As')
    vasp_calc.set_computer(code.get_computer())
    vasp_calc.set_queue_name(queue)
    vasp_calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine': 20})
    vasp_calc.label = 'Test VASP run'
    vasp_calc.store_all()
    vasp_calc.submit()


@run_example.command()
@example_param_set
def simple(paw_family, import_from, queue, code, computer):
    load_dbenv_if_not_loaded()
    from aiida.orm import CalculationFactory, Code
    if import_from:
        import_paws(import_from, paw_family)
    paw_si, = get_paws(paw_family, 'Si')

    vasp_calc = CalculationFactory('vasp.vasp')()
    vasp_calc.use_structure(create_structure_Si())
    vasp_calc.use_kpoints(create_kpoints())
    vasp_calc.use_parameters(create_params_simple())
    code = Code.get_from_string('{}@{}'.format(code, computer))
    vasp_calc.use_code(code)
    vasp_calc.use_paw(paw_si, 'Si')
    vasp_calc.set_computer(code.get_computer())
    vasp_calc.set_queue_name(queue)
    vasp_calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine': 20})
    vasp_calc.label = 'Test VASP run'
    vasp_calc.store_all()
    vasp_calc.submit()


def load_dbenv_if_not_loaded():
    from aiida import load_dbenv, is_dbenv_loaded
    if not is_dbenv_loaded():
        load_dbenv()


def get_data_cls(descriptor):
    load_dbenv_if_not_loaded()
    from aiida.orm import DataFactory
    return DataFactory(descriptor)


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
    return param_cls(dict={
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
    return param_cls(dict={'prec': 'NORMAL', 'encut': 200, 'ediff': 1e-8, 'ialgo': 38, 'ismear': 0.0, 'sigma': 0.1})


def import_paws(folder_path, family_name):
    load_dbenv_if_not_loaded()
    from aiida.orm import DataFactory
    paw_cls = DataFactory('vasp.paw')
    paw_cls.import_family(folder_path, familyname=family_name, family_desc='Test family')


def get_paws(family_name, *symbols):
    load_dbenv_if_not_loaded()
    from aiida.orm import DataFactory
    paw_cls = DataFactory('vasp.paw')
    try:
        paws = [paw_cls.load_paw(family=family_name, symbol=symbol)[0] for symbol in symbols]
        return paws
    except ValueError as err:
        click.echo(repr(err), err=True)
        raise ValueError('give a valid family or import a new one (run with --help)')


if __name__ == '__main__':
    run_example()
