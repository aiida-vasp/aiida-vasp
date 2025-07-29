import pytest
from aiida import orm
from ase.build import bulk

from aiida_vasp.workchains.v2.common import builder_updater as bup


def test_vasp_builder_updater(aiida_profile_clean, vasp_code, upload_potcar, potcar_family_name):
    structure = orm.StructureData(ase=bulk('MgO', 'rocksalt', 5.0)).store()
    vasp_code.store()
    upd = bup.VaspBuilderUpdater()
    upd.apply_preset(structure, code='vasp@localhost')
    assert upd.builder.structure == structure
    assert upd.builder.parameters['incar']['algo'] == 'normal'
    assert upd.builder.options['resources']['tot_num_mpiprocs'] == 1
    assert upd.builder.code == vasp_code
    assert upd.builder.settings.get_dict() == {}
    assert upd.builder.kpoints_spacing.value == 0.05
    assert upd.builder.potential_mapping.get_dict() == {'Mg': 'Mg_pv', 'O': 'O'}
    # The default family is not uploaded - this should remain None
    assert upd.builder.potential_family is None

    upd.set_potential_family(potcar_family_name)
    assert upd.builder.potential_family.value == potcar_family_name


def test_vasp_relax_updater(aiida_profile_clean, vasp_code):
    structure = orm.StructureData(ase=bulk('MgO', 'rocksalt', 5.0)).store()
    vasp_code.store()
    upd = bup.VaspRelaxUpdater()
    upd.apply_preset(structure, code='vasp@localhost')
    assert upd.builder.structure == structure
    assert upd.builder.vasp.parameters['incar']['algo'] == 'normal'
    assert upd.builder.vasp.options['resources']['tot_num_mpiprocs'] == 1
    assert upd.builder.vasp.code == vasp_code
    assert upd.builder.vasp.settings.get_dict() == {}
    assert upd.builder.vasp.kpoints_spacing.value == 0.05
    assert upd.builder.vasp.potential_family is None
    assert upd.builder.vasp.potential_mapping.get_dict() == {'Mg': 'Mg_pv', 'O': 'O'}
    assert upd.builder.relax_settings.get_dict()['algo'] == 'cg'


@pytest.mark.parametrize('hybrid', [True, False])
def test_vasp_band_updater(aiida_profile_clean, vasp_code, hybrid):
    structure = orm.StructureData(ase=bulk('MgO', 'rocksalt', 5.0)).store()
    vasp_code.store()
    if hybrid:
        upd = bup.VaspHybridBandUpdater()
    else:
        upd = bup.VaspBandUpdater()
    upd.apply_preset(structure, code='vasp@localhost')
    assert upd.builder.structure == structure
    assert upd.builder.scf.parameters['incar']['algo'] == 'normal'
    assert upd.builder.scf.options['resources']['tot_num_mpiprocs'] == 1
    assert upd.builder.scf.code == vasp_code
    assert upd.builder.scf.settings.get_dict() == {}
    assert upd.builder.scf.kpoints_spacing.value == 0.05
    assert upd.builder.scf.potential_family is None
    assert upd.builder.scf.potential_mapping.get_dict() == {'Mg': 'Mg_pv', 'O': 'O'}

    if hybrid:
        upd = bup.VaspHybridBandUpdater()
    else:
        upd = bup.VaspBandUpdater()
    upd.apply_preset(structure, code='vasp@localhost', run_relax=True)

    assert upd.builder.structure == structure
    assert upd.builder.relax.vasp.parameters['incar']['algo'] == 'normal'
    assert upd.builder.relax.vasp.options['resources']['tot_num_mpiprocs'] == 1
    assert upd.builder.relax.vasp.code == vasp_code
    assert upd.builder.relax.vasp.settings.get_dict() == {}
    assert upd.builder.relax.vasp.kpoints_spacing.value == 0.05
    assert upd.builder.relax.vasp.potential_family is None
    assert upd.builder.relax.vasp.potential_mapping.get_dict() == {'Mg': 'Mg_pv', 'O': 'O'}


def test_vasp_neb_updater(aiida_profile_clean, vasp_code):
    structure = orm.StructureData(ase=bulk('MgO', 'rocksalt', 5.0)).store()
    vasp_code.store()
    upd = bup.VaspNEBUpdater()
    upd.apply_preset(structure, structure, code='vasp@localhost')
