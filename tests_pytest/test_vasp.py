from gaas_sample import get_gaas_process_inputs


def test_vasp(configure_with_daemon, get_gaas_process_inputs, assert_finished):
    from aiida.work.run import run

    process, inputs = get_gaas_process_inputs(calculation_string='vasp.vasp')
    output, pid = run(
        process,
        _return_pid=True,
        **inputs
    )
    assert_finished(pid)
