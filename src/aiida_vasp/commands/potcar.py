"""
Commands for the potential interface.

-------------------------------------
Commandline util for dealing with potcar files.
"""

# pylint: disable=import-outside-toplevel
import click
import tabulate
from aiida.cmdline.utils.decorators import with_dbenv
from click_spinner import spinner as cli_spinner

from aiida_vasp.commands import options
from aiida_vasp.utils.aiida_utils import cmp_load_verdi_data, get_data_class

VERDI_DATA = cmp_load_verdi_data()


@VERDI_DATA.group('vasp-potcar')
def potcar():
    """Top level command for handling VASP POTCAR files."""


def try_grab_description(ctx, param, value):
    """
    Try to get the description from an existing group if it's not given.

    This is a click parameter callback.
    """
    potcar_data_cls = get_data_class('vasp.potcar')
    group_name = ctx.params['name']
    existing_groups = potcar_data_cls.get_potcar_groups()
    existing_group_names = [group.label for group in existing_groups]
    if not value:
        if group_name in existing_group_names:
            return potcar_data_cls.get_potcar_group(group_name).description
        raise click.MissingParameter('A new group must be given a description.', param=param)
    return value


def detect_old_style_groups():
    """Check for the existence of old style groups and prompt the user"""
    from aiida.orm import Group, QueryBuilder

    from aiida_vasp.data.potcar import OLD_POTCAR_FAMILY_TYPE, PotcarGroup

    qdb = QueryBuilder()
    qdb.append(Group, filters={'type_string': OLD_POTCAR_FAMILY_TYPE}, project=['label'])
    all_old_groups = [qres[0] for qres in qdb.all()]
    not_migrated = []
    for group_label in all_old_groups:
        qdb = QueryBuilder()
        qdb.append(PotcarGroup, filters={'label': {'==': group_label}})
        count = qdb.count()
        if count == 0:
            not_migrated.append(group_label)
    if any(not_migrated):
        click.echo(
            (
                'Some of the old style POTCAR family groups are not migrated. '
                "Please run command 'verdi data vasp-potcar migratefamilies.\n",
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

    potcar_data_cls = get_data_class('vasp.potcar')
    with cli_spinner():
        num_found, num_added, num_uploaded = potcar_data_cls.upload_potcar_family(
            path, name, description, stop_if_existing=stop_if_existing, dry_run=dry_run
        )

    click.echo(f'POTCAR files found: {num_found}. New files uploaded: {num_uploaded}, Added to Family: {num_added}')
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

    potcar_data_cls = get_data_class('vasp.potcar')
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
    potcar_data_cls = get_data_class('vasp.potcar')

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
