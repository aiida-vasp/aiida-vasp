"""
Test the export command.
"""

from pathlib import Path

import pytest
from click.testing import CliRunner

from aiida_vasp.commands.tools import tools


def run_cmd(command=None, args=None, **kwargs):
    """Run verdi data vasp.tools <command> [args]."""
    runner = CliRunner()
    params = args or []
    if command:
        params.insert(0, command)
    return runner.invoke(tools, params, **kwargs)


# TODO - add tests for other commands
# Combine export test with workflow execution test to save time
@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_uploadfamily_withpath(aiida_profile_clean, tmp_path, run_vasp_process):
    """
    Test export vasp calculation
    """

    _, node = run_vasp_process(test_case='exit_codes/converged')
    result = run_cmd(
        'export',
        args=[str(node.pk), str(Path(tmp_path) / str(node.pk))],
    )
    assert result.exit_code == 0

    assert (Path(tmp_path) / f'{node.pk}/INCAR').is_file()
    assert (Path(tmp_path) / f'{node.pk}/OUTCAR').is_file()
    assert (Path(tmp_path) / f'{node.pk}/vasprun.xml').is_file()
