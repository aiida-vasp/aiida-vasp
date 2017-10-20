"""Vasp calculation integration test"""
from gaas_sample import get_gaas_process_inputs  # pylint: disable=unused-import


def test_vasp(configure_with_daemon, get_gaas_process_inputs, assert_finished):  # pylint: disable=unused-argument,redefined-outer-name
    """Vasp calculation integration test"""

    from aiida.work.run import run
    from aiida.orm import DataFactory

    process, inputs = get_gaas_process_inputs(
        calculation_string='vasp.vasp', parameters={
            'ncore': 4
        })
    inputs.settings = DataFactory('parameter')(dict={
        'ADDITIONAL_RETRIEVE_LIST': ['CHGCAR']
    })
    output, pid = run(process, _return_pid=True, **inputs)
    assert_finished(pid)
    assert all(filename in output['retrieved'].get_folder_list()
               for filename in
               ['OUTCAR', 'EIGENVAL', 'DOSCAR', 'CHGCAR', 'vasprun.xml'])
