from gaas_sample import get_gaas_process_inputs


def test_vasp2w90(configure_with_daemon, sample, get_gaas_process_inputs, assert_finished):
    from aiida.work.run import run
    from aiida.orm import DataFactory

    process, inputs = get_gaas_process_inputs(
        calculation_string='vasp.vasp2w90'
    )
    charge_density = DataFactory('vasp.chargedensity')()
    charge_density.set_file(sample('GaAs/CHGCAR'))
    inputs.charge_density = charge_density
    output, pid = run(
        process,
        _return_pid=True,
        **inputs
    )
    print(output)
    assert_finished(pid)
