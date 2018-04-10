"""Unittests for PoscarIo"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.poscar import PoscarParser


@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_poscar_io(fresh_aiida_env, vasp_structure_poscar):
    poscario = vasp_structure_poscar
    assert poscario.potentials_order == ['In', 'As', 'In_d', 'As']


@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_parse_poscar(fresh_aiida_env, vasp_structure):
    """Parse a reference POSCAR file with the PoscarParser and compare the result to a reference structure."""
    path = data_path('poscar', 'POSCAR')
    parser = PoscarParser(path, 'POSCAR', None)
    result = parser.get_quantity(None, 'structure', {})
    structure = vasp_structure
    assert result['structure'].cell == structure.cell
    assert result['structure'].get_site_kindnames() == structure.get_site_kindnames()
    assert result['structure'].sites[2].position == structure.sites[2].position
