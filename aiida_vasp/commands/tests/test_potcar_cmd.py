"""Unit tests for vasp-potcar command family."""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
from __future__ import absolute_import
from __future__ import print_function
import os
import pytest

from py import path as py_path  # pylint: disable=no-member,no-name-in-module
from click.testing import CliRunner
from monty.collections import AttrDict
from packaging import version

from aiida_vasp.commands.potcar import potcar
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.utils.fixtures.environment import fresh_aiida_env
from aiida_vasp.utils.fixtures.data import potcar_family, POTCAR_FAMILY_NAME, temp_pot_folder


@pytest.fixture
def cmd_params(temp_pot_folder):
    """Common building blocks for ``uploadfamily`` calls."""
    params = AttrDict()
    params.POTCAR_PATH = temp_pot_folder.strpath
    params.FAMILY_NAME = POTCAR_FAMILY_NAME
    params.PATH_OPTION = '--path={}'.format(params.POTCAR_PATH)
    params.NAME_OPTION = '--name={}'.format(params.FAMILY_NAME)
    params.DESC_OPTION = '--description="This is a test POTCAR family"'
    return params


def run_cmd(command=None, args=None, **kwargs):
    """Run verdi data vasp-potcar <command> [args]."""
    runner = CliRunner()
    params = args or []
    if command:
        params.insert(0, command)
    return runner.invoke(potcar, params, **kwargs)


def test_no_subcmd():
    result = run_cmd()
    assert not result.exception


def test_uploadfamily_withpath(fresh_aiida_env, cmd_params):
    """Upload the test potcar family and check it is there."""

    result = run_cmd('uploadfamily', [cmd_params.PATH_OPTION, cmd_params.NAME_OPTION, cmd_params.DESC_OPTION])

    potcar_cls = get_data_class('vasp.potcar')

    assert not result.exception
    assert potcar_cls.exists(element='In')
    assert potcar_cls.exists(element='Ga')
    assert [g.label for g in potcar_cls.get_potcar_groups()] == [cmd_params.FAMILY_NAME]


def test_uploadfamily_tar(fresh_aiida_env, cmd_params):
    """Give a tar file as the source."""
    path_option = '--path={}'.format(py_path.local(cmd_params.POTCAR_PATH).join('Ga.tar'))
    result = run_cmd('uploadfamily', [path_option, cmd_params.NAME_OPTION, cmd_params.DESC_OPTION])
    potcar_cls = get_data_class('vasp.potcar')

    assert not result.exception
    assert potcar_cls.exists(element='Ga')
    assert [g.label for g in potcar_cls.get_potcar_groups()] == [cmd_params.FAMILY_NAME]


def test_uploadfamily_inworkdir(fresh_aiida_env, cmd_params):
    """Upload the test potcar family from the working env."""

    potcar_dir = py_path.local(cmd_params.POTCAR_PATH)
    old_work_dir = potcar_dir.chdir()

    result = run_cmd('uploadfamily', [cmd_params.NAME_OPTION, cmd_params.DESC_OPTION])

    potcar_cls = get_data_class('vasp.potcar')

    assert not result.exception
    assert potcar_cls.exists(element='In')
    assert [g.label for g in potcar_cls.get_potcar_groups()] == [cmd_params.FAMILY_NAME]

    old_work_dir.chdir()


def test_uploadfamily_again(fresh_aiida_env, potcar_family, cmd_params):
    """
    Re-upload a potcar family.

    Test:
        * Does not require description
        * Must succeed
        * Adds no nodes
        * Adds no groups
    """
    from aiida.orm import Node, Group
    from aiida.orm.querybuilder import QueryBuilder

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


def test_uploadfamily_dryrun(fresh_aiida_env, cmd_params):
    """Make sure --dry-run does not affect the db."""
    from aiida.orm import Node, Group
    from aiida.orm.querybuilder import QueryBuilder

    node_qb = QueryBuilder(path=[Node])
    node_count = node_qb.count()
    group_qb = QueryBuilder(path=[Group])
    group_count = group_qb.count()

    result = run_cmd('uploadfamily', [cmd_params.PATH_OPTION, cmd_params.NAME_OPTION, cmd_params.DESC_OPTION, '--dry-run'])

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


def test_listfamilies_nofilter(fresh_aiida_env, potcar_family):
    """Test typical usecases without filtering."""
    result = run_cmd('listfamilies')
    assert not result.exception
    assert potcar_family in result.output

    family_group = get_data_class('vasp.potcar').get_potcar_group(potcar_family)
    result = run_cmd('listfamilies', ['--description'])
    assert not result.exception
    assert 'Description' in result.output
    assert family_group.description in result.output


def test_listfamilies_filtering(fresh_aiida_env, potcar_family):
    """Test filtering families by elements & symbols."""
    result = run_cmd('listfamilies', ['--element', 'In', '--element', 'As'])
    assert potcar_family in result.output

    result = run_cmd('listfamilies', ['--element', 'In', '--element', 'U235'])
    assert potcar_family not in result.output

    result = run_cmd('listfamilies', ['--symbol', 'In_d'])
    assert potcar_family in result.output

    result = run_cmd('listfamilies', ['--symbol', 'In_d', '--symbol', 'In_s'])
    assert potcar_family not in result.output

    result = run_cmd('listfamilies', ['--symbol', 'In_d', '--element', 'As'])
    assert potcar_family in result.output

    result = run_cmd('listfamilies', ['--symbol', 'In_d', '--element', 'U235'])
    assert potcar_family not in result.output


def test_exportfamilies(fresh_aiida_env, potcar_family, tmpdir):
    """Test exporting potcar family."""
    result = run_cmd('exportfamily', ['--name', potcar_family, '--path', str(tmpdir)])
    assert not result.exception
    export_path = tmpdir.join(potcar_family)
    assert export_path.isdir()
    assert export_path.exists()

    new_dir = tmpdir.join('new_dir')
    result = run_cmd('exportfamily', ['--dry-run', '--name', potcar_family, '--path', str(new_dir)])
    assert not result.exception
    assert not new_dir.exists()

    result = run_cmd('exportfamily', ['--as-archive', '--name', potcar_family, '--path', str(tmpdir)])
    assert not result.exception
    export_path = tmpdir.join(potcar_family + '.tar.gz')
    assert export_path.isfile()
    assert export_path.exists()

    new_arch = tmpdir.join('export.tar.gz')
    result = run_cmd('exportfamily', ['--dry-run', '--as-archive', '--name', potcar_family, '--path', str(new_arch)])
    assert not result.exception
    assert not new_arch.exists()


def test_call_from_vasp():
    """Test if the verdi potcar data command works."""

    import subprocess
    output = subprocess.check_output(['verdi', 'data', 'vasp-potcar', '--help'], universal_newlines=True)
    assert 'Usage: verdi data vasp-potcar' in output  # pylint: disable=unsupported-membership-test
