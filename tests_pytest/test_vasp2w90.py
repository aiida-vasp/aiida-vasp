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
    from aiida.orm.data.base import List

    process, inputs = get_gaas_process_inputs(
        calculation_string='vasp.vasp2w90',
        parameters={
            'ncore': 1,
            'isym': -1,
            'icharg': 11,
            'lwave': False
        })
    charge_density = DataFactory('vasp.chargedensity')()
    charge_density.set_file(sample('GaAs/CHGCAR'))
    inputs.charge_density = charge_density

    wannier_parameters = DataFactory('parameter')(dict=dict(
        dis_num_iter=1000,
        num_bands=24,
        num_iter=0,
        num_wann=14,
        spinors=True,
    ))
    inputs.wannier_parameters = wannier_parameters

    wannier_projections = List()
    wannier_projections.extend(['Ga : s; px; py; pz', 'As : px; py; pz'])
    inputs.wannier_projections = wannier_projections

    output, pid = run(process, _return_pid=True, **inputs)
    assert all(key in output
               for key in [
                   'retrieved', 'kpoints', 'wannier_parameters',
                   'wannier_kpoints', 'wannier_projections'
               ])
    assert_finished(pid)
