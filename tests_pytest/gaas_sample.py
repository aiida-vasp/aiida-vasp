"""GaAs calculation integration test fixture"""
try:
    from collections import ChainMap
except ImportError:
    from chainmap import ChainMap

import pytest
from ase.io.vasp import read_vasp


@pytest.fixture
# pylint: disable=protected-access
def get_gaas_process_inputs(configure, get_process_inputs, sample):  # pylint: disable=unused-argument
    """Test fixture: sets up a calculation to be tested"""

    def inner(calculation_string, parameters=None):  # pylint: disable=missing-docstring
        from aiida.orm import DataFactory

        parameters = parameters if parameters else {}

        process, inputs = get_process_inputs(calculation_string=calculation_string, code_string='vasp')

        structure = DataFactory('structure')()
        structure.set_ase(read_vasp(sample('GaAs/POSCAR')))
        inputs.structure = structure

        paw_cls = DataFactory('vasp.paw')
        ga_paw = paw_cls.load_paw(family='pbe', symbol='Ga')[0]
        as_paw = paw_cls.load_paw(family='pbe', symbol='As')[0]
        inputs.paw = {'Ga': ga_paw, 'As': as_paw}

        kpoints = DataFactory('array.kpoints')()
        kpoints.set_kpoints_mesh([8, 8, 8])
        inputs.kpoints = kpoints

        parameters_input = DataFactory('parameter')(
            dict=ChainMap(parameters, dict(
                nbands=24,
                gga='PE',
                gga_compat=False,
                ismear=0,
                lorbit=11,
                lsorbit=True,
                sigma=0.05,
            )))
        inputs.parameters = parameters_input
        inputs._options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 12}
        inputs._options.withmpi = True
        inputs._options.queue_name = 'dphys_compute'
        inputs._options.max_wallclock_seconds = 600

        return process, inputs

    return inner
