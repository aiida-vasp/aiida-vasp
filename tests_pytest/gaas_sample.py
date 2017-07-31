import pytest
from ase.io.vasp import read_vasp

@pytest.fixture
def get_gaas_process_inputs(configure, get_process_inputs, sample):
    def inner(calculation_string):
        from aiida.orm import DataFactory

        process, inputs = get_process_inputs(
            calculation_string=calculation_string,
            code_string='vasp'
        )

        structure = DataFactory('structure')()
        structure.set_ase(read_vasp(sample('GaAs/POSCAR')))
        inputs.structure = structure

        Paw = DataFactory('vasp.paw')
        Ga_paw = Paw.load_paw(family='pbe', symbol='Ga')[0]
        As_paw = Paw.load_paw(family='pbe', symbol='As')[0]
        inputs.paw = {'Ga': Ga_paw, 'As': As_paw}

        kpoints = DataFactory('array.kpoints')()
        kpoints.set_kpoints_mesh([8, 8, 8])
        inputs.kpoints = kpoints

        parameters = DataFactory('parameter')(dict=dict(
            npar=8,
            gga='PE',
            gga_compat=False,
            ismear=0,
            lorbit=11,
            lsorbit=True,
            sigma=0.05,
        ))
        inputs.parameters = parameters
        inputs._options.resources = {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 16
        }
        inputs._options.withmpi = True
        inputs._options.queue_name = 'dphys_compute'

        return process, inputs
    return inner
