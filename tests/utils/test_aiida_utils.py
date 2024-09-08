"""Test aiida_utils functionss."""

import pytest

from aiida_vasp.utils.aiida_utils import (
    get_current_user,
)
from aiida_vasp.utils.pmg import PymatgenAdapator, get_incar, get_kpoints, get_outcar, get_vasprun


def test_get_current_user(fresh_aiida_env):
    """Assert that get_current_user returns a user."""
    user = get_current_user()
    assert user.pk
    assert user.first_name == ''
    assert user.last_name == ''
    assert user.email


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_pmg_adaptor(fresh_aiida_env, tmp_path, run_vasp_process):
    """
    Test export vasp calculation
    """

    _, node = run_vasp_process()

    adapt = PymatgenAdapator(node)
    assert adapt.vasprun

    assert 'vasprun' in adapt.pmg_objects

    with PymatgenAdapator(node) as adapt:
        adapt.vasprun

    assert 'vasprun_dict' in adapt.cache
    assert 'pmg_cache' in node.base.extras.all

    assert get_incar(node)
    assert get_kpoints(node)
    assert get_vasprun(node)
    assert get_outcar(node)
