"""Vasp2W90 calculation integration test"""
from __future__ import print_function

from gaas_sample import get_gaas_process_inputs  # pylint: disable=unused-import


def test_vasp2w90(
        configure_with_daemon,  # pylint: disable=unused-argument
        sample,
        get_gaas_process_inputs,  # pylint: disable=redefined-outer-name
        assert_finished):
    """Vasp2W90 unit test"""
    from aiida.work.run import run
    from aiida.orm import DataFactory

    process, inputs = get_gaas_process_inputs(
        calculation_string='vasp.vasp2w90')
    charge_density = DataFactory('vasp.chargedensity')()
    charge_density.set_file(sample('GaAs/CHGCAR'))
    inputs.charge_density = charge_density
    output, pid = run(process, _return_pid=True, **inputs)
    assert all(key in output
               for key in [
                   'retrieved',
                   'kpoints',
                   'wannier_parameters',
                   'wannier_kpoints',
               ])
    assert_finished(pid)
