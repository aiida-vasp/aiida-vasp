"""
Commands for the potential interface.

Commandline util for dealing with potcar files.
"""

# pylint: disable=import-outside-toplevel
from pathlib import Path

import click
import tabulate
from aiida.cmdline.utils.decorators import with_dbenv
from click_spinner import spinner as cli_spinner

from aiida_vasp.commands import cmd_aiida_vasp, options

# from aiida_vasp.utils.aiida_utils import cmp_load_verdi_data

# VERDI_DATA = cmp_load_verdi_data()

FUNCTIONAL_CHOICES = [
    'PBE',
    'PBE_52',
    'PBE_52_W_HASH',
    'PBE_54',
    'PBE_54_W_HASH',
    'PBE_64',
    'LDA',
    'LDA_52',
    'LDA_52_W_HASH',
    'LDA_54',
    'LDA_54_W_HASH',
    'LDA_64',
    'PW91',
    'LDA_US',
    'PW91_US',
    'Perdew_Zunger81',
]


# @VERDI_DATA.group('vasp.potcar')
@cmd_aiida_vasp.group('potcar')
def potcar() -> None:
    """Top level command for handling VASP POTCAR files."""


def try_grab_description(ctx, param, value):
    """
    Try to get the description from an existing group if it's not given.

    This is a click parameter callback.
    """
    from aiida_vasp.data.potcar import PotcarData

    potcar_data_cls = PotcarData
    group_name = ctx.params['name']
    existing_groups = potcar_data_cls.get_potcar_groups()
    existing_group_names = [group.label for group in existing_groups]
    if not value:
        if group_name in existing_group_names:
            return potcar_data_cls.get_potcar_group(group_name).description
        raise click.MissingParameter('A new group must be given a description.', param=param)
    return value


def detect_old_style_groups() -> None:
    """Check for the existence of old style groups and prompt the user"""

    from aiida import orm

    from aiida_vasp.data.potcar import OLD_POTCAR_FAMILY_TYPE, PotcarGroup

    qdb = orm.QueryBuilder()
    qdb.append(orm.Group, filters={'type_string': OLD_POTCAR_FAMILY_TYPE}, project=['label'])
    all_old_groups = [qres[0] for qres in qdb.all()]
    not_migrated = []
    for group_label in all_old_groups:
        qdb = orm.QueryBuilder()
        qdb.append(PotcarGroup, filters={'label': {'==': group_label}})
        count = qdb.count()
        if count == 0:
            not_migrated.append(group_label)
    if any(not_migrated):
        click.echo(
            (
                'Some of the old style POTCAR family groups are not migrated. '
                "Please run command 'aiida-vasp potcar migratefamilies.\n",
                f'The missing groups are: {not_migrated}.',
            )
        )


@potcar.command()
@options.PATH(
    help='Path to a folder or archive containing the POTCAR files. '
    'You can supply the archive that you downloaded from the VASP server. '
    'The path does not need to be specified, if that is the case, the current path is used.'
)
@options.FAMILY_NAME()
@options.DESCRIPTION(help='A description for the family.', callback=try_grab_description)
@click.option(
    '--stop-if-existing', is_flag=True, help='An option to abort when encountering a previously uploaded POTCAR file.'
)
@options.DRY_RUN()
@with_dbenv()
def uploadfamily(path, name, description, stop_if_existing, dry_run):
    """Upload a family of VASP potcar files."""

    from aiida_vasp.data.potcar import PotcarData

    potcar_data_cls = PotcarData
    with cli_spinner():
        num_found, num_added, num_uploaded = potcar_data_cls.upload_potcar_family(
            path, name, description, stop_if_existing=stop_if_existing, dry_run=dry_run
        )

    click.echo(f'POTCAR files found: {num_found}. New files uploaded: {num_uploaded}, Added to Family: {num_added}')
    if dry_run:
        click.echo('No files were uploaded due to --dry-run.')


@potcar.command()
@click.argument('family_name')
def listsymbols(family_name):
    """List available symbols in a POTCAR family group."""
    from aiida import orm
    from aiida.cmdline.utils import echo

    from aiida_vasp.data.potcar import PotcarGroup

    group: PotcarGroup = orm.load_group(family_name)
    symbols = [(node.symbol, node.original_file_name) for node in group.nodes]
    click.echo(f"Symbols in family '{family_name}':")
    for symbol, fpath in symbols:
        click.echo(f'- {symbol:<20} -> {fpath}')

    duplicated, _ = group.get_duplicated_symbols()
    if duplicated:
        echo.echo_warning(f'Duplicated symbols found in group {family_name}: {duplicated}')


@potcar.command()
@click.argument('family_name')
@click.option('--show-each', is_flag=True, help='Show the resolution of each symbol')
def integrity(family_name, show_each):
    """Check the integrity of a POTCAR family"""
    from aiida import orm
    from aiida.cmdline.utils import echo

    from aiida_vasp.data.potcar import PotcarGroup

    group: PotcarGroup = orm.load_group(family_name)
    click.echo(f'Group: {family_name} with {group.count()} potcars')
    duplicated, resolved = group.get_duplicated_symbols()
    if duplicated:
        echo.echo_warning(f'Duplicated symbols found in group {family_name}: {duplicated}')
        for symbol, fname in resolved.items():
            echo.echo(f'{symbol} -> {fname}')
    else:
        echo.echo_success('No duplicated symbols found')
    matched = group.get_matched_set()
    if matched:
        echo.echo_success(f'This group matches an known set `{matched}`')
    else:
        echo.echo('This group does not match any known set')
    if show_each:
        echo.echo('Match of individual POTCARs:')
        identities = group.get_potcar_identity()
        for key, value in identities.items():
            echo.echo(f'{key} ->  {value}')


@potcar.command()
@click.argument('family_name', required=False)
@click.option('--dryrun', required=False, is_flag=True)
def fix_inconsistent_symbols(family_name, dryrun):
    """Fix inconsistent families"""
    from uuid import uuid4

    from aiida import orm
    from aiida.cmdline.utils import echo

    from aiida_vasp.data.potcar import PotcarData, PotcarGroup, UniquenessError

    if family_name is not None:
        groups = [orm.load_group(family_name)]
    else:
        groups = PotcarGroup.collection.all()
        groups = [g for g in groups if ('.transient-' not in g.label) and (not g.label.endswith('.backup'))]
    echo.echo(f'Groups to check: {[group.label for group in groups]}')
    for group in groups:
        duplicated, _ = group.get_duplicated_symbols()
        symbols = [node.symbol for node in group.nodes]
        problematic_nodes = {}
        with cli_spinner():
            for node in group.nodes:
                info = node.check_and_fix_inconsistent_potcar_symbol()
                if info is not None:
                    problematic_nodes[node.pk] = info
        # Now check if the result matches with the duplicated info
        for dup_symbol, ndup in duplicated.items():
            node_pks = [key for key, value in problematic_nodes.items() if value['stored_symbol'] == dup_symbol]
            if len(node_pks) != ndup - 1:
                echo.echo(
                    f'Inconsistent POTCARs for symbol {dup_symbol} - expect {ndup - 1} inconsistent nodes but '
                    f'{len(node_pks)} found'
                )
                for node_pk in node_pks:
                    echo.echo(
                        f' - {node_pk} symbol: {problematic_nodes[node_pk]["stored_symbol"]} '
                        f'original_file_name: {problematic_nodes[node_pk]["original_filename"]}'
                    )
                # echo.echo_critical('Abort')
        # Fix - first let's collect the updated nodes
        updated_nodes = []
        for pk in problematic_nodes.keys():
            node: PotcarData = orm.load_node(pk)
            info = node.check_and_fix_inconsistent_potcar_symbol(fix=True)
            new_node = info['updated_node']
            updated_nodes.append(new_node)

            # Should not happen, but worth checking
            if new_node.symbol in symbols:
                echo.echo_critical(
                    f'Symbol `{new_node.symbol}` of the new node {new_node} already exists in the group {group}!'
                    ' Something is seriously wrong.'
                )

        if not updated_nodes:
            echo.echo(f'No inconsistent nodes found in group {group.label}')
            continue
        echo.echo('Nodes to be removed:')
        for pk in problematic_nodes.keys():
            node = orm.load_node(pk)
            echo.echo(f' - {node} INCONSISTENT symbol: {node.symbol} original_file_name: {node.original_file_name} ')
        echo.echo('\n\n')
        echo.echo('Nodes to be added:')
        for node in updated_nodes:
            echo.echo(f' - {node} UPDATED symbol: {node.symbol} original_file_name: {node.original_file_name} ')
        if dryrun:
            echo.echo('Dryrun - no changes made')

        if click.confirm(f'Proceed to fix group {group.label}?', default=False) and not dryrun:
            # Create a backup group
            new_group = PotcarGroup(label=group.label + f'.transient-{str(uuid4())[:8]}')
            new_group.store()
            new_group.add_nodes(list(group.nodes))
            echo.echo_info(f'Temporary new Group {new_group} created')

            # Store the updated nodes
            for node in updated_nodes:
                if not node.is_stored:
                    try:
                        node.store()
                    except UniquenessError as e:
                        echo.echo_error(f'Failed to store node {node}: {e}')
            echo.echo_info('New PotcarData nodes stored.')

            # Remove the problematic nodes from the original group
            new_group.remove_nodes([orm.load_node(pk) for pk in problematic_nodes.keys()])
            echo.echo_info('Inconsistent PotcarData nodes removed.')

            # Add the updated nodes to the original group
            new_group.add_nodes(updated_nodes)
            echo.echo_info('Updated PotcarData nodes added.')

            # Check again
            dup, _ = new_group.get_duplicated_symbols()
            if len(dup) > 0:
                echo.echo(f'Operation failed - there are still duplicated symbols in the group {new_group}')
                echo.echo_critical(f'Abort. {group} remain unchanged')

            # Swap the label
            old_nodes = list(group.nodes)
            backup_group = PotcarGroup(label=group.label + '.backup')
            backup_group.store()
            backup_group.add_nodes(old_nodes)
            group.remove_nodes(old_nodes)
            group.add_nodes(list(new_group.nodes))

            echo.echo_success(f'Successfully fixed group {group.label}')
            echo.echo(f'You can now delete transient group {new_group} with `verdi group delete {new_group.pk}`.')


@potcar.command()
@click.option(
    '--functional',
    help='Name of the functional to be used for the POTCAR files.',
    type=click.Choice(FUNCTIONAL_CHOICES),
    default='PBE',
)
@options.FAMILY_NAME()
@options.DESCRIPTION(help='A description for the family.', callback=try_grab_description)
@click.option(
    '--stop-if-existing', is_flag=True, help='An option to abort when encountering a previously uploaded POTCAR file.'
)
@options.DRY_RUN()
@with_dbenv()
def upload_from_pymatgen(functional, name, description, stop_if_existing, dry_run):
    """
    Upload a family of VASP potcar files from pymatgen

    If you have pymatgen installed and configured to locate the correct VASP POTCAR files,
    you can use this command to upload a family of VASP potcar files into aiida-vasp.
    """
    from pymatgen.io.vasp.inputs import SETTINGS, PotcarSingle

    from aiida_vasp.data.potcar import PotcarData
    from aiida_vasp.utils.pmg import convert_pymatgen_potcar_folder, temporary_folder

    funcdir = PotcarSingle.functional_dir[functional]
    pmg_vasp_psp_dir = SETTINGS.get('PMG_VASP_PSP_DIR')
    if pmg_vasp_psp_dir is None:
        raise click.Abort(
            'PMG_VASP_PSP_DIR is not set, please set it in your .pmgrc.yaml file or set the environment variable'
        )
    source_folder = f'{pmg_vasp_psp_dir}/{funcdir}'
    if not Path(source_folder).exists():
        raise click.Abort(f'The source folder {source_folder} does not exist.')
    with cli_spinner():
        # Convert the pymatgen potcar folder to a temporary folder with the same structure as the VASP potcar folder
        with temporary_folder() as temp_folder:
            convert_pymatgen_potcar_folder(source_folder, temp_folder)
            # Try to upload from this folder
            num_found, num_added, num_uploaded = PotcarData.upload_potcar_family(
                temp_folder, name, description, stop_if_existing=stop_if_existing, dry_run=dry_run
            )
            click.echo(
                f'POTCAR files found: {num_found}. New files uploaded: {num_uploaded}, Added to Family: {num_added}'
            )
            if dry_run:
                click.echo('No files were uploaded due to --dry-run.')


@potcar.command()
@click.option(
    '-e', '--element', multiple=True, help='Filter for families containing potentials for all given elements.'
)
@click.option('-s', '--symbol', multiple=True, help='Filter for families containing potentials for all given symbols.')
@click.option('-d', '--description', is_flag=True, help='Also show the description.')
@with_dbenv()
def listfamilies(element, symbol, description):
    """List available families of VASP potcar files."""
    detect_old_style_groups()

    from aiida_vasp.data.potcar import PotcarData

    potcar_data_cls = PotcarData
    groups = potcar_data_cls.get_potcar_groups(filter_elements=element, filter_symbols=symbol)

    table = [['Family', 'Num Potentials']]
    if description:
        table[0].append('Description')
    for group in groups:
        row = [group.label, len(group.nodes)]
        if description:
            row.append(group.description)
        table.append(row)
    if len(table) > 1:
        click.echo(tabulate.tabulate(table, headers='firstrow'))
        click.echo()
    elif element or symbol:
        click.echo('No POTCAR family contains all given elements and symbols.')
    else:
        click.echo('No POTCAR family available.')


@potcar.command()
@options.PATH(type=click.Path(exists=False), help='Path to location of the exported POTCAR family.')
@options.FAMILY_NAME()
@options.DRY_RUN(help='Only display what would be exported.')
@click.option('-z', '--as-archive', is_flag=True, help='Create a compressed archive (.tar.gz) instead of a folder.')
@click.option('-v', '--verbose', is_flag=True, help='Print the names of all created files.')
@with_dbenv()
def exportfamily(path, name, dry_run, as_archive, verbose):
    """Export a POTCAR family into a compressed tar archive or folder."""
    from aiida_vasp.data.potcar import PotcarData

    potcar_data_cls = PotcarData

    if not as_archive:
        files = potcar_data_cls.export_family_folder(name, path, dry_run)
        if verbose:
            click.echo(tabulate.tabulate([[i] for i in files], headers=['Files written:']))
    else:
        archive, files = potcar_data_cls.export_family_archive(name, path, dry_run)
        if verbose:
            click.echo(tabulate.tabulate([[i] for i in files], headers=[f'Files added to archive {archive}:']))

    click.echo(f'{len(files)} POTCAR files exported.')
    if dry_run:
        click.echo('Nothing written due to "--dry-run"')


@potcar.command()
@with_dbenv()
def migratefamilies():
    """
    Migrate the type_string associated with the potcar family groups.

    Previously, these groups has type_string: data.vasp.potcar.family.
    Since AiiDA 1.2, groups used by plugins should be defined by subclass and entrypoint names.
    This commands recreates the old style group using the ``PotcarGroup`` class.
    """

    from aiida_vasp.data.potcar import migrate_potcar_group

    migrate_potcar_group()
