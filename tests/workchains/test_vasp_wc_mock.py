"""
WorkChain test using mock code

The tests here uses mock-vasp to simulate running VASP.
To generate the test data, set the following environment variables:

- `MOCK_VASP_POTCAR_PATH`: path to the directory containing the POTCAR files
- `MOCK_VASP_VASP_CMD`: command to run VASP

When generating the mock data, make sure to add:

```
export MOCK_VASP_UPLOAD_PREFIX=<test_case_name>
```

to the `custom_scheduler_commands` so the uploaded folder has tag prefix to it.
Otherwise, it will be impossible to link the test data to the test cases.

The test data are stored in the `<root>/tests/test_data` folder.

Make sure you unset the environmental variables and rerun the tests to check it works as before.
"""

from aiida import orm


def test_silicon_sp(mock_potcars, mock_vasp_strict, builder_updater):
    """Test running a VASP workchain on silicon using the mock code."""
    from ase.build import bulk

    si = bulk('Si', 'diamond', 5.4)
    si_node = orm.StructureData(ase=si)

    upd = builder_updater

    upd.apply_preset(si_node, code='mock-vasp@localhost')
    upd.set_options(custom_scheduler_commands='export MOCK_VASP_UPLOAD_PREFIX=mock_silicon_sp')
    results = upd.run_get_node()
    # Add prefix to the registry folder
    assert results.node.is_finished_ok


def test_silicon_relax(mock_potcars, mock_vasp_strict, builder_updater):
    """Test running a VASP workchain on silicon using the mock code."""
    from ase.build import bulk

    si = bulk('Si', 'diamond', 5.4)
    si_node = orm.StructureData(ase=si)

    from aiida_vasp.workchains import VaspRelaxWorkChain

    upd = VaspRelaxWorkChain.get_builder_updater(code='mock-vasp@localhost')
    upd.apply_preset(si_node)
    upd.set_options(custom_scheduler_commands='export MOCK_VASP_UPLOAD_PREFIX=mock_silicon_relax')
    results = upd.run_get_node()
    # Add prefix to the registry folder
    assert results.node.is_finished_ok


def test_silicon_converge(mock_potcars, mock_vasp_strict):
    """Test running a VASP workchain on silicon using the mock code."""
    from ase.build import bulk

    si = bulk('Si', 'diamond', 5.4)
    si_node = orm.StructureData(ase=si)

    from aiida_vasp.workchains import VaspConvergenceWorkChain

    upd = VaspConvergenceWorkChain.get_builder_updater(code='mock-vasp@localhost')
    upd.apply_preset(si_node)
    upd.set_conv_settings(cutoff_stop=400.0, kspacing_stop=0.06)
    # Add prefix to the registry folder
    upd.set_options(custom_scheduler_commands='export MOCK_VASP_UPLOAD_PREFIX=mock_silicon_convergence')
    results = upd.run_get_node()
    assert results.node.is_finished_ok


def test_silicon_band(mock_potcars, mock_vasp_strict):
    """Test running a VASP workchain on silicon using the mock code."""
    from ase.build import bulk

    si = bulk('Si', 'diamond', 5.4)
    si_node = orm.StructureData(ase=si)

    from aiida_vasp.workchains import VaspBandsWorkChain

    upd = VaspBandsWorkChain.get_builder_updater(code='mock-vasp@localhost')
    upd.apply_preset(si_node)
    # Add prefix to the registry folder
    upd.set_options(custom_scheduler_commands='export MOCK_VASP_UPLOAD_PREFIX=mock_silicon_bands')
    results = upd.run_get_node()
    assert results.node.is_finished_ok


def test_silicon_band_hybrid(mock_potcars, mock_vasp_strict):
    """Test the hybrid (split-path) SCF  band structure workchain"""

    from ase.build import bulk

    from aiida_vasp.workchains import VaspHybridBandsWorkChain

    si = bulk('Si', 'diamond', 5.4)
    si_node = orm.StructureData(ase=si)
    upd = VaspHybridBandsWorkChain.get_builder_updater(code='mock-vasp@localhost')
    upd.apply_preset(si_node, run_relax=True)
    upd.set_band_settings(kpoints_per_split=120)
    # Add prefix to the registry folder
    upd.set_options(custom_scheduler_commands='export MOCK_VASP_UPLOAD_PREFIX=mock_silicon_hybrid')
    results = upd.run_get_node()
    assert results.node.is_finished_ok
