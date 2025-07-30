"""
Unit tests for vasp.potcar command family.
"""

import os
import subprocess
from pathlib import Path

import pytest
from aiida.common import AttributeDict
from aiida.orm import Group, Node
from aiida.orm.querybuilder import QueryBuilder
from click.testing import CliRunner

from aiida_vasp.commands.potcar import potcar
from aiida_vasp.data.potcar import PotcarData, PotcarGroup


@pytest.fixture
def cmd_params(temp_pot_folder, potcar_family_name):
    """Common building blocks for ``uploadfamily`` calls."""
    params = AttributeDict()
    params.POTCAR_PATH = str(temp_pot_folder)
    params.FAMILY_NAME = potcar_family_name
    params.PATH_OPTION = f'--path={params.POTCAR_PATH}'
    params.NAME_OPTION = f'--name={params.FAMILY_NAME}'
    params.DESC_OPTION = '--description="This is a test POTCAR family"'
    return params


def run_cmd(command=None, args=None, **kwargs):
    """Run verdi data vasp.potcar <command> [args]."""
    runner = CliRunner()
    params = args or []
    if command:
        params.insert(0, command)
    return runner.invoke(potcar, params, **kwargs)


def test_no_subcmd():
    result = run_cmd('--help')
    assert result.exception is None


def test_uploadfamily_withpath(aiida_profile_clean, cmd_params):
    """Upload the test potcar family and check it is there."""

    result = run_cmd(
        'uploadfamily',
        [cmd_params.PATH_OPTION, cmd_params.NAME_OPTION, cmd_params.DESC_OPTION],
    )

    potcar_cls = PotcarData

    assert not result.exception
    assert potcar_cls.exists(element='In')
    assert potcar_cls.exists(element='Ga')
    assert [g.label for g in potcar_cls.get_potcar_groups()] == [cmd_params.FAMILY_NAME]


def test_uploadfamily_tar(aiida_profile_clean, cmd_params):
    """Give a tar file as the source."""
    path_option = f'--path={Path(cmd_params.POTCAR_PATH) / "Ga.tar"!s}'
    result = run_cmd('uploadfamily', [path_option, cmd_params.NAME_OPTION, cmd_params.DESC_OPTION])
    potcar_cls = PotcarData

    assert not result.exception
    assert potcar_cls.exists(element='Ga')
    assert [g.label for g in potcar_cls.get_potcar_groups()] == [cmd_params.FAMILY_NAME]


def test_uploadfamily_inworkdir(aiida_profile_clean, cmd_params):
    """Upload the test potcar family from the working env."""

    potcar_dir = Path(cmd_params.POTCAR_PATH)
    old_work_dir = Path().cwd()
    os.chdir(str(potcar_dir))
    assert str(potcar_dir) == str(Path().cwd())

    result = run_cmd('uploadfamily', [cmd_params.NAME_OPTION, cmd_params.DESC_OPTION])

    potcar_cls = PotcarData

    assert not result.exception
    assert potcar_cls.exists(element='In')
    assert [g.label for g in potcar_cls.get_potcar_groups()] == [cmd_params.FAMILY_NAME]

    os.chdir(str(old_work_dir))

    assert str(old_work_dir) == str(Path().cwd())


def test_uploadfamily_again(aiida_profile_clean, upload_potcar, cmd_params):
    """
    Re-upload a potcar family.

    Test:
        * Does not require description
        * Must succeed
        * Adds no nodes
        * Adds no groups
    """

    node_qb = QueryBuilder(path=[Node])
    node_count = node_qb.count()
    group_qb = QueryBuilder(path=[Group])
    group_count = group_qb.count()

    result = run_cmd('uploadfamily', [cmd_params.PATH_OPTION, cmd_params.NAME_OPTION])

    assert not result.exception

    node_qb = QueryBuilder(path=[Node])
    assert node_count == node_qb.count()
    group_qb = QueryBuilder(path=[Group])
    assert group_count == group_qb.count()


def test_uploadfamily_dryrun(aiida_profile_clean, cmd_params):
    """Make sure --dry-run does not affect the db."""
    node_qb = QueryBuilder(path=[Node])
    node_count = node_qb.count()
    group_qb = QueryBuilder(path=[Group])
    group_count = group_qb.count()

    result = run_cmd(
        'uploadfamily',
        [
            cmd_params.PATH_OPTION,
            cmd_params.NAME_OPTION,
            cmd_params.DESC_OPTION,
            '--dry-run',
        ],
    )

    assert not result.exception

    node_qb = QueryBuilder(path=[Node])
    assert node_count == node_qb.count()
    group_qb = QueryBuilder(path=[Group])
    assert group_count == group_qb.count()


def test_listfamilies_existence():
    """Make sure the listfamilies subcommand exists."""
    result = run_cmd('listfamilies')
    assert not result.exception
    assert result.output


def test_listfamilies_nofilter(aiida_profile_clean, upload_potcar, potcar_family_name):
    """Test typical usecases without filtering."""
    result = run_cmd('listfamilies')
    assert not result.exception
    assert potcar_family_name in result.output

    family_group = PotcarData.get_potcar_group(potcar_family_name)
    result = run_cmd('listfamilies', ['--description'])
    assert not result.exception
    assert 'Description' in result.output
    assert family_group.description in result.output


def test_listfamilies_filtering(aiida_profile_clean, upload_potcar, potcar_family_name):
    """Test filtering families by elements & symbols."""
    result = run_cmd('listfamilies', ['--element', 'In', '--element', 'As'])
    assert potcar_family_name in result.output

    result = run_cmd('listfamilies', ['--element', 'In', '--element', 'U235'])
    assert potcar_family_name not in result.output

    result = run_cmd('listfamilies', ['--symbol', 'In_d'])
    assert potcar_family_name in result.output

    result = run_cmd('listfamilies', ['--symbol', 'In_d', '--symbol', 'In_s'])
    assert potcar_family_name not in result.output

    result = run_cmd('listfamilies', ['--symbol', 'In_d', '--element', 'As'])
    assert potcar_family_name in result.output

    result = run_cmd('listfamilies', ['--symbol', 'In_d', '--element', 'U235'])
    assert potcar_family_name not in result.output


def test_exportfamilies(aiida_profile_clean, upload_potcar, potcar_family_name, tmp_path):
    """Test exporting potcar family."""
    result = run_cmd('exportfamily', ['--name', potcar_family_name, '--path', str(tmp_path)])
    assert not result.exception
    export_path = tmp_path / potcar_family_name
    assert export_path.is_dir()
    assert export_path.exists()

    new_dir = tmp_path / 'new_dir'
    result = run_cmd('exportfamily', ['--dry-run', '--name', potcar_family_name, '--path', new_dir])
    assert not result.exception
    assert not new_dir.exists()

    result = run_cmd('exportfamily', ['--as-archive', '--name', potcar_family_name, '--path', tmp_path])
    assert not result.exception
    name = potcar_family_name + '.tar.gz'
    export_path = tmp_path / name
    assert export_path.is_file()
    assert export_path.exists()

    new_arch = tmp_path / 'export.tar.gz'
    result = run_cmd(
        'exportfamily',
        ['--dry-run', '--as-archive', '--name', potcar_family_name, '--path', new_arch],
    )
    assert not result.exception
    assert not new_arch.exists()


def test_call_from_vasp():
    """Test if the verdi potcar data command works."""

    output = subprocess.check_output(['verdi', 'data', 'vasp.potcar', '--help'], universal_newlines=True)
    assert 'Usage: verdi data vasp.potcar' in output  # pylint: disable=unsupported-membership-test


def test_migrate_command(aiida_profile_clean, legacy_potcar_family):
    """Test the migration command"""

    legacy_name, legacy_group_class = legacy_potcar_family
    legacy_group = legacy_group_class.collection.get(label=legacy_name)
    run_cmd('migratefamilies')
    migrated = PotcarGroup.collection.get(label=legacy_name)

    uuids_original = {node.uuid for node in legacy_group.nodes}
    uuids_migrated = {node.uuid for node in migrated.nodes}
    assert uuids_migrated == uuids_original
