"""Unittests for PoscarIo"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
import pytest

from aiida_vasp.utils.fixtures import *


@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_poscar_io(fresh_aiida_env, vasp_structure_poscar):
    poscario = vasp_structure_poscar
    assert poscario.potentials_order == ['In', 'As', 'In_d', 'As']
