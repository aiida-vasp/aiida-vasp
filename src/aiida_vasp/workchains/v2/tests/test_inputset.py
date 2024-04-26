"""
Test for input set specifications
"""

import pytest
from aiida_vasp.workchains.v2.inputset.base import InputSet
from aiida_vasp.workchains.v2.inputset.vaspsets import VASPInputSet
from ase.build import bulk

# pylint:disable=redefined-outer-name,unused-argument


@pytest.fixture
def fe_atoms():
    """Get a Fe atoms"""
    return bulk('Fe', 'fcc', 5.0)


@pytest.fixture
def mgo_atoms():
    """Get a MgO atoms"""
    return bulk('MgO', 'rocksalt', 5.0)


def test_base(fe_atoms):
    """Base test case"""
    iset = InputSet('MITRelaxSet', fe_atoms, overrides={'ediff': 1.0, 'nsw': None})

    out = iset.get_input_dict()
    assert out['ediff'] == 1.0
    assert out['ibrion'] == 2
    assert 'nsw' not in out


def test_vasp(fe_atoms):
    """Test VASP inputsets"""
    iset = VASPInputSet('MITRelaxSet', fe_atoms, overrides={'ediff': 1.0, 'nsw': None, 'ldautype': 3})

    out = iset.get_input_dict()
    assert out['ediff'] == 1.0
    assert out['ibrion'] == 2
    assert out['magmom'] == [5]
    assert out['ldauu'] == [4.0]
    assert out['ldauj'] == [0.0]
    assert out['ldaul'] == [2]
    assert out['ldautype'] == 3
    assert out['ldau'] is True
    assert 'nsw' not in out


def test_kpoints(aiida_profile, fe_atoms):
    """Test generating kpoints"""
    inset = VASPInputSet('MITRelaxSet', fe_atoms)
    kpoints = inset.get_kpoints(0.05)
    assert kpoints.get_kpoints_mesh()[0][0] == 7
