"""Vasp2W90 calculation integration test"""
from __future__ import print_function

from gaas_sample import get_gaas_process_inputs  # pylint: disable=unused-import


def test_vasp2w90(sample, get_gaas_process_inputs, assert_finished):  # pylint: disable=redefined-outer-name
    """Vasp2W90 unit test"""
    from aiida.work.run import run
    from aiida.orm import DataFactory

    process, inputs = get_gaas_process_inputs(
        calculation_string='vasp.vasp2w90',
        parameters={'ncore': 1,
                    'isym': -1,
                    'icharg': 11,
                    'lwave': False})
    charge_density = DataFactory('vasp.chargedensity')()
    charge_density.set_file(sample('GaAs/CHGCAR'))
    inputs.charge_density = charge_density
    output, pid = run(process, _return_pid=True, **inputs)
    print(output)
    assert_finished(pid)
