"""Unit tests for vasp-potcar command family"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name
import os
from py import path as py_path  # pylint: disable=no-member,no-name-in-module
from click.testing import CliRunner

from aiida_vasp.commands.potcar import potcar
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.utils.fixtures.environment import aiida_env, fresh_aiida_env
from aiida_vasp.utils.fixtures.data import potcar_family, POTCAR_FAMILY_NAME

POTCAR_PATH = data_path('potcar')
FAMILY_NAME = POTCAR_FAMILY_NAME
PATH_OPTION = '--path={}'.format(POTCAR_PATH)
NAME_OPTION = '--name={}'.format(FAMILY_NAME)
DESC_OPTION = '--description="This is a test POTCAR family"'


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


def test_uploadfamily_withpath(fresh_aiida_env):
    """Upload the test potcar family and check it is there."""

    result = run_cmd('uploadfamily', [PATH_OPTION, NAME_OPTION, DESC_OPTION])

    potcar_cls = get_data_class('vasp.potcar')

    assert not result.exception
    assert potcar_cls.exists(element='In')
    assert potcar_cls.exists(element='Ga')
    assert [g.name for g in potcar_cls.get_potcar_groups()] == [FAMILY_NAME]


def test_uploadfamily_tar(fresh_aiida_env):
    """Give a tar file as the source"""
    path_option = '--path={}'.format(py_path.local(POTCAR_PATH).join('Ga.tar'))
    result = run_cmd('uploadfamily', [path_option, NAME_OPTION, DESC_OPTION])
    potcar_cls = get_data_class('vasp.potcar')

    print result.output

    assert not result.exception
    assert potcar_cls.exists(element='Ga')
    assert [g.name for g in potcar_cls.get_potcar_groups()] == [FAMILY_NAME]


def test_uploadfamily_inworkdir(fresh_aiida_env):
    """Upload the test potcar family from the working env."""

    potcar_dir = py_path.local(POTCAR_PATH)
    old_work_dir = potcar_dir.chdir()

    result = run_cmd('uploadfamily', [NAME_OPTION, DESC_OPTION])

    potcar_cls = get_data_class('vasp.potcar')

    assert not result.exception
    assert potcar_cls.exists(element='In')
    assert [g.name for g in potcar_cls.get_potcar_groups()] == [FAMILY_NAME]

    old_work_dir.chdir()


def test_uploadfamily_again(potcar_family):
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

    result = run_cmd('uploadfamily', [PATH_OPTION, NAME_OPTION])

    assert not result.exception

    node_qb = QueryBuilder(path=[Node])
    assert node_count == node_qb.count()
    group_qb = QueryBuilder(path=[Group])
    assert group_count == group_qb.count()


def test_uploadfamily_dryrun(fresh_aiida_env):
    """Make sure --dry-run does not affect the db"""
    from aiida.orm import Node, Group
    from aiida.orm.querybuilder import QueryBuilder

    node_qb = QueryBuilder(path=[Node])
    node_count = node_qb.count()
    group_qb = QueryBuilder(path=[Group])
    group_count = group_qb.count()

    result = run_cmd('uploadfamily', [PATH_OPTION, NAME_OPTION, DESC_OPTION, '--dry-run'])

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


def test_listfamilies_nofilter(potcar_family):
    """Test typical usecases without filtering."""
    result = run_cmd('listfamilies')
    assert not result.exception
    assert potcar_family in result.output

    family_group = get_data_class('vasp.potcar').get_potcar_group(potcar_family)
    result = run_cmd('listfamilies', ['--with-description'])
    assert not result.exception
    assert 'Description' in result.output
    assert family_group.description in result.output


def test_listfamilies_filtering(potcar_family):
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


def test_exportfamilies(potcar_family, tmpdir):
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
