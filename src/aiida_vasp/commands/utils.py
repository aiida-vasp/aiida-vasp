"""
Auxsiliary functions for AiiDA VASP vasp.launch command line tools.
"""

from pathlib import Path

import click
from aiida import orm
from aiida.common.exceptions import NotExistent
from aiida.plugins import DataFactory
from ase.io import read

from aiida_vasp.commands.option_parser import process_dict_option
from aiida_vasp.common.builder_updater import VaspBuilderUpdater
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser


def load_structure(structure_path: str | Path) -> orm.StructureData:
    """Load a structure from various file formats."""
    structure_path = Path(structure_path)
    if not structure_path.exists():
        try:
            structure = orm.load_node(structure_path)
        except NotExistent:
            raise click.ClickException(f'Structure file not found: {structure_path} nor it is a valid node identifier')
        # Try a node
        return structure

    # Try to determine format from extension
    extension = structure_path.suffix.lower()

    if extension in ['.cif']:
        # Load CIF file
        cif_data = DataFactory('core.cif')
        cif_node, _ = cif_data.get_or_create(str(structure_path.absolute()))
        # Convert to StructureData
        structure = orm.StructureData(ase=cif_node.get_ase())
    elif extension in ['.poscar', '.contcar', '.vasp'] or structure_path.name.upper() in ['POSCAR', 'CONTCAR']:
        # Load POSCAR/CONTCAR file
        with open(structure_path, 'r', encoding='utf8') as handler:
            poscar_parser = PoscarParser(handler=handler)
        poscar_dict = poscar_parser.structure
        structure = orm.StructureData()
        structure.set_cell(poscar_dict['unitcell'])
        for site in poscar_dict['sites']:
            structure.append_atom(position=site['position'], symbols=site['symbol'], name=site['kind_name'])
    else:
        # Try to use ASE for other formats (XYZ, etc.)
        try:
            atoms = read(str(structure_path))
            structure = orm.StructureData(ase=atoms)
        except Exception as e:
            raise click.ClickException(f'Could not load structure from {structure_path}: {e}')

    # Set a label based on filename
    structure.label = structure_path.stem
    return structure


def setup_calculation_options(options, resources, max_wallclock_seconds, num_machines, tot_num_mpiprocs):
    """Setup computational resources from various options."""
    options_dict = {}
    if options:
        options_dict.update(process_dict_option(options))
    if resources:
        options_dict['resources'] = process_dict_option(resources)
    if max_wallclock_seconds:
        options_dict['max_wallclock_seconds'] = max_wallclock_seconds
    if 'resources' not in options_dict:
        options_dict['resources'] = {}
    if num_machines:
        options_dict['resources']['num_machines'] = num_machines
    if tot_num_mpiprocs:
        options_dict['resources']['tot_num_mpiprocs'] = tot_num_mpiprocs
    return options_dict


def apply_additional_updates(upd: VaspBuilderUpdater, additional_overrides: dict):
    """
    Apply additional overrides to the builder updater by using the set_xxx methods.
    """
    if not additional_overrides:
        return

    click.echo(f'Loading input overrides from: {additional_overrides}')

    # Apply other overrides
    for key, value in additional_overrides.items():
        if hasattr(upd, f'set_{key}'):
            method = getattr(upd, f'set_{key}')
            if isinstance(value, dict):
                method(**value)
            else:
                method(value)
            click.echo(f'Applied {key} overrides')


def handle_calculation_submission(
    upd: VaspBuilderUpdater, run_directly: bool, group: str, alias: str | None = None
) -> orm.ProcessNode:
    """Handle calculation submission and group assignment."""
    # Submit or run the calculation
    if not run_directly:
        click.echo('Submitting calculation to daemon...')
        result = upd.submit()
        click.echo(f'Submitted calculation with PK: {result.pk}')
        click.echo(f'UUID: {result.uuid}')
        process_node = result
    else:
        click.echo('Running calculation directly...')
        result = upd.run_get_node()
        click.echo(f'Calculation completed with PK: {result.node.pk}')
        click.echo(f'UUID: {result.node.uuid}')
        if result.node.is_finished_ok:
            click.echo('Calculation finished successfully!')
        else:
            click.echo(f'Calculation failed with exit status: {result.node.exit_status}')
            if result.node.exit_message:
                click.echo(f'Exit message: {result.node.exit_message}')
        process_node = result.node

    # Add to group if specified
    if group:
        try:
            calc_group, created = orm.Group.collection.get_or_create(label=group)
            # Using GroupPathX - we give the resulting node an alias if specified
            if alias:
                try:
                    from aiida_grouppathx.pathx import GroupPathX
                except ImportError:
                    raise ImportError('aiida-grouppathx is required for setting an alias. Please install it first.')
                GroupPathX(group).add_node(process_node, alias=alias)
                alias_action = f' with alias "{alias}"'
            else:
                # Just add the node normally
                calc_group.add_nodes([process_node])
                alias_action = ''
            action = 'Created' if created else 'Added to existing'
            click.echo(f"{action} group '{group}' and added calculation{alias_action}")
        except Exception as e:
            click.echo(f"Warning: Could not add to group '{group}': {e}")

    return process_node
