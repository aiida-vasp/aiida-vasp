import click
import numpy


@click.command()
def test_vasp():
    load_dbenv_if_not_loaded()
    from aiida.orm import CalculationFactory, Code
    vasp_calc = CalculationFactory('vasp.vasp')()
    vasp_calc.use_structure(create_structure())
    vasp_calc.use_kpoints(create_kpoints())
    vasp_calc.use_parameters(create_params())
    code = Code.get_from_string('vasp@monch')
    vasp_calc.use_code(code)
    vasp_calc.set_computer(code.get_computer())
    vasp_calc.set_resources({
        'num_machines': 1,
        'num_mpiprocs_per_machine': 20
    })
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


def create_structure():
    structure_cls = get_data_cls('structure')
    structure = structure_cls(
        cell=numpy.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * 6.058, )
    structure.append_atom(position=(0, 0, 0), symbols='In')
    structure.append_atom(position=(0.25, 0.25, 0.25), symbols='As')
    return structure


def create_kpoints():
    kpoints_cls = get_data_cls('array.kpoints')
    return kpoints_cls(kpoints_mesh=[8, 8, 8])


def create_params():
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


if __name__ == '__main__':
    test_vasp()
