"""Commandline util for dealing with potcar files"""
import click
import tabulate

from aiida.cmdline.commands import data_cmd

from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.commands import options


@data_cmd.group('vasp-potcar')
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
    existing_group_names = [group.name for group in existing_groups]
    if not value:
        if group_name in existing_group_names:
            return potcar_data_cls.get_potcar_group(group_name).description
        else:
            raise click.MissingParameter('A new group must be given a description.', param=param)
    return value


@potcar.command()
@options.PATH(help='Path to a folder or archive containing the POTCAR files')
@options.FAMILY_NAME()
@options.DESCRIPTION(help='A description for the family', callback=try_grab_description)
@click.option('--stop-if-existing', is_flag=True, help='Abort when encountering a previously uploaded POTCAR file')
@options.DRY_RUN()
def uploadfamily(path, name, description, stop_if_existing, dry_run):
    """Upload a family of VASP potcar files."""

    potcar_data_cls = get_data_class('vasp.potcar')
    num_found, num_added, num_uploaded = potcar_data_cls.upload_potcar_family(
        path, name, description, stop_if_existing=stop_if_existing, dry_run=dry_run)

    click.echo('POTCAR files found: {}. New files uploaded: {}, Added to Family: {}'.format(num_found, num_uploaded, num_added))
    if dry_run:
        click.echo('No files were uploaded due to --dry-run.')


@potcar.command()
@click.option('-e', '--element', multiple=True, help='Filter for families containing potentials for all given elements.')
@click.option('-s', '--symbol', multiple=True, help='Filter for families containing potentials for all given symbols.')
@click.option('-d', '--with-description', is_flag=True)
def listfamilies(element, symbol, with_description):
    """List available families of VASP potcar files."""

    potcar_data_cls = get_data_class('vasp.potcar')
    groups = potcar_data_cls.get_potcar_groups(filter_elements=element, filter_symbols=symbol)

    table = [['Family', 'Num Potentials']]
    if with_description:
        table[0].append('Description')
    for group in groups:
        row = [group.name, len(group.nodes)]
        if with_description:
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
            click.echo(tabulate.tabulate([[i] for i in files], headers=['Files added to archive {}:'.format(archive)]))

    click.echo('{} POTCAR files exported.'.format(len(files)))
    if dry_run:
        click.echo('Nothing written due to "--dry-run"')
