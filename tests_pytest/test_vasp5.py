import sys
from ase.io.vasp import read_vasp

def test_vasp5(configure_with_daemon, sample, get_process_inputs, assert_finished):
    from aiida.orm import Code, DataFactory, CalculationFactory
    from aiida.work.run import run
    from aiida_vasp.calcs.maker import VaspMaker

    process, inputs = get_process_inputs(
        calculation_string='vasp.vasp5',
        code_string='vasp'
    )

    structure = DataFactory('structure')()
    structure.set_ase(read_vasp(sample('GaAs/POSCAR')))
    inputs.structure = structure

    Paw = DataFactory('vasp.paw')
    Ga_paw = Paw.load_paw(family='pbe', symbol='Ga')[0]
    print(type(Ga_paw))
    inputs.paw = {
        'Ga': Ga_paw,
        'As': Paw.load_paw(family='pbe', symbol='As')[0]
    }

    kpoints = DataFactory('array.kpoints')()
    kpoints.set_kpoints_mesh([8, 8, 8])
    inputs.kpoints = kpoints

    settings = DataFactory('parameter')(dict=dict(
        npar=8,
        gga='PE',
        gga_compat=False,
        ismear=0,
        lorbit=11,
        lsorbit=True,
        sigma=0.05,
    ))
    inputs.settings = settings
    inputs._options.resources = {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 16
    }
    inputs._options.withmpi = True
    inputs._options.queue_name = 'dphys_compute'
    output, pid = run(
        process,
        _return_pid=True,
        **inputs
    )
    assert_finished(pid)
