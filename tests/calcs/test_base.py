"""Unit tests for aiida_vasp.calcs.base."""

# pylint: disable=unused-import,unused-argument,redefined-outer-name, import-outside-toplevel
import os

import pytest
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.engine.utils import instantiate_process
from aiida.manage.manager import get_manager

from aiida_vasp.calcs.base import VaspCalcBase


def test_generate_base_calc(base_calc):
    """Test that it is possible to start the generation of an instance of the base calculation class."""

    manager = get_manager()
    runner = manager.get_runner()
    inputs = AttributeDict()
    metadata = AttributeDict({'options': {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}})
    inputs.metadata = metadata

    # Need more input, so we will get an AttributeError because we do not pass code
    with pytest.raises(AttributeError):
        instantiate_process(runner, VaspCalcBase, **inputs)

    # Let us try again with a code attached, thus should pass
    inputs.code = base_calc.inputs.code
    instantiate_process(runner, VaspCalcBase, **inputs)


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('cif', 'mesh')], indirect=True)
def test_remote_copy_restart_folder_without_settings(vasp_calc, vasp_inputs, vasp_code, tmp_path):
    """A restart folder without a ``settings`` input must not raise.

    ``vasp.relax`` with ``reuse=True`` attaches its previous remote folder to
    the final single-point as ``restart_folder``; protocols (and users) do not
    have to provide ``settings``, and this used to except at submit time with
    ``AttributeError: 'AttributesFrozendict' object has no attribute
    'settings'``.
    """
    (tmp_path / 'CHGCAR').write_text('fake-chgcar', encoding='utf8')
    (tmp_path / 'WAVECAR').write_text('fake-wavecar', encoding='utf8')
    inputs = vasp_inputs(settings=None)
    assert 'settings' not in inputs
    inputs.restart_folder = orm.RemoteData(computer=vasp_code.computer, remote_path=str(tmp_path)).store()
    calc = vasp_calc(inputs=inputs)
    copy_list = calc.remote_copy_restart_folder()
    assert sorted(os.path.basename(entry[1]) for entry in copy_list) == ['CHGCAR', 'WAVECAR']
