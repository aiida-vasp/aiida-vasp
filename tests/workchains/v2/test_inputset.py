"""
Test for input set specifications
"""

import numpy as np
import pytest
from aiida import orm
from ase.build import bulk

from aiida_vasp.workchains.v2.inputset.base import InputSet
from aiida_vasp.workchains.v2.inputset.pmgset import PymatgenInputSet
from aiida_vasp.workchains.v2.inputset.vaspsets import VASPInputSet

try:
    # Import SETTINGS
    from pymatgen.core import SETTINGS

    # Check equivalence with pymatgen
    from pymatgen.io.vasp.sets import MPRelaxSet
except ImportError:
    MPRelaxSet = None

# pylint:disable=redefined-outer-name,unused-argument


@pytest.fixture
def fe_atoms(aiida_profile):
    """Get a Fe atoms"""
    return orm.StructureData(ase=bulk('Fe', 'fcc', 5.0))


@pytest.fixture
def feo_atoms(aiida_profile):
    """Get a FeO atoms in rocksalt structure"""
    return orm.StructureData(ase=bulk('FeO', 'rocksalt', 5.0))


@pytest.fixture
def mgo_atoms(aiida_profile):
    """Get a MgO atoms"""
    return orm.StructureData(ase=bulk('MgO', 'rocksalt', 5.0))


def test_base(fe_atoms):
    """Base test case"""
    iset = InputSet('MITRelaxSet', overrides={'ediff': 1.0, 'nsw': None})

    out = iset.get_input_dict(fe_atoms)
    assert out['ediff'] == 1.0
    assert out['ibrion'] == 2
    assert 'nsw' not in out


def test_vasp(fe_atoms):
    """Test VASP inputsets"""
    iset = VASPInputSet('MITRelaxSet', overrides={'ediff': 1.0, 'nsw': None, 'ldautype': 3})

    out = iset.get_input_dict(fe_atoms)
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
    inset = VASPInputSet('MITRelaxSet')
    kpoints = inset.get_kpoints(fe_atoms, 0.05)
    assert kpoints.get_kpoints_mesh()[0][0] == 7


@pytest.mark.skipif(MPRelaxSet is None, reason='pymatgen is not installed')
def test_pmg_kpoints(aiida_profile, fe_atoms):
    """Test generating kpoints"""
    inset = PymatgenInputSet('MPRelaxSet')
    kpoints = inset.get_kpoints(fe_atoms)
    assert kpoints.get_kpoints_mesh()[0][0] == 7
    assert kpoints.get_kpoints_mesh()[1][0] == 0

    inset = PymatgenInputSet('MPRelaxSet')
    scaled_fe_atoms = fe_atoms.get_ase()
    scaled_fe_atoms.set_cell(scaled_fe_atoms.get_cell() * 1.2)
    # Breaks the symmetry so it is not face centred
    scaled_fe_atoms.cell += np.random.rand(3, 3) * 0.1
    scaled_fe_atoms = orm.StructureData(ase=scaled_fe_atoms)
    kpoints = inset.get_kpoints(scaled_fe_atoms)
    assert kpoints.get_kpoints_mesh()[0][0] == 6
    assert kpoints.get_kpoints_mesh()[1][0] == -0.5

    inset = PymatgenInputSet('MP24RelaxSet')
    kpoints = inset.get_kpoints(fe_atoms)
    assert kpoints is None
    kspacing = inset.get_kpoints_spacing(fe_atoms)
    assert kspacing == 0.22 / np.pi / 2


@pytest.mark.skipif(MPRelaxSet is None, reason='pymatgen is not installed')
def test_pmg(fe_atoms, feo_atoms):
    """Test pymatgen inputsets"""
    vpmgset = PymatgenInputSet('MPRelaxSet', overrides={'ediff': 1.0, 'nsw': None, 'ldautype': 3})
    out = vpmgset.get_input_dict(fe_atoms)
    assert 'nsw' not in out
    assert out['ediff'] == 1.0
    assert out['ldautype'] == 3
    assert 'ldau' not in out
    assert out['ismear'] == -5

    vpmgset = PymatgenInputSet('MPRelaxSet', overrides={'ediff': 1.0, 'nsw': None})
    out = vpmgset.get_input_dict(feo_atoms)
    assert 'ldau' in out

    psp_dir = SETTINGS.get('PMG_VASP_PSP_DIR')
    for structure in feo_atoms, fe_atoms:
        pmgset = MPRelaxSet(structure.get_pymatgen())
        vpmgset = PymatgenInputSet('MPRelaxSet')
        ref_dict = dict(pmgset.incar)
        ref_dict.pop('icharg', None)
        ref_dict.pop('istart', None)
        ref_dict.pop('kspacing', None)
        out = vpmgset.get_input_dict(structure)
        assert ref_dict == {key.upper(): value for key, value in out.items()}
        assert pmgset.kpoints.kpts[0] == tuple(vpmgset.get_kpoints(structure).get_kpoints_mesh()[0])
        if psp_dir:
            assert {p.element: p.symbol for p in pmgset.potcar} == vpmgset.get_pp_mapping(structure)
