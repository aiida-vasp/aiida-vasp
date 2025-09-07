"""
Tool for importing existing calculations into the database
"""

import os
import tempfile
from pathlib import Path

import click

from . import cmd_aiida_vasp


@cmd_aiida_vasp.command('import')
@click.argument('path')
@click.option('--code', '-c', required=False, help='Code to use for the calculation.')
@click.option('--potential-family', '-f', required=False, help='Potential family to use if POTCAR is missing.')
@click.option('--potential-mapping', '-m', required=False, help='JSON string mapping elements to potentials.')
@click.option('--include-wavecar', is_flag=True, help='Include WAVECAR file if present.')
@click.option('--include-chgcar', is_flag=True, help='Include CHGCAR file if present.')
@click.option('--stdout-file', default='vasp_output', help='Name of the stdout file to look for.')
@click.option('--label', '-l', required=False, help='Label for the imported calculation.')
@click.option('--description', '-d', required=False, help='Description for the imported calculation.')
@click.option('--group', '-g', required=False, help='Group to add the calculation to.')
@click.option('--submit-daemon', is_flag=True, help='Submit to daemon instead of running directly.', default=False)
@click.option('--quiet', '-q', is_flag=True, help='Suppress output except for errors.')
@click.option('--yes', '-y', is_flag=True, help='Automatically confirm all prompts.')
def import_calc(
    code,
    path,
    potential_family,
    potential_mapping,
    include_wavecar,
    include_chgcar,
    stdout_file,
    label,
    description,
    group,
    submit_daemon,
    quiet,
    yes,
):
    """
    Import an existing calculation into the database.
    The calculation will be imported as a `VaspCalculation` with inputs/outputs connected.
    TODO: automatically link to existing StructureData node using hashing mechanism.
    """
    import json

    from aiida import orm
    from aiida.common.exceptions import NotExistent

    from aiida_vasp.calcs.immigrant import VaspCalcImporter

    if not quiet:
        click.echo(f'Importing calculation from: {path}')

    # Check the existence of the calculation folder
    calc_path = Path(path)
    is_remote = False

    if not calc_path.exists():
        # This might be a remote path, we'll handle it later
        is_remote = True
        if not quiet:
            click.echo(f'Local path {path} does not exist - treating as remote path')
    else:
        # Check for required VASP input files
        required_files = ['INCAR', 'POSCAR', 'KPOINTS']
        missing_files = [f for f in required_files if not (calc_path / f).exists()]

        if missing_files:
            click.echo(f'Error: Missing required files: {missing_files}', err=True)
            raise click.Abort()

        # Check for POTCAR
        has_potcar = (calc_path / 'POTCAR').exists()
        if not has_potcar and not potential_family:
            click.echo('Warning: POTCAR not found and no potential family specified.', err=True)
            click.echo('You must provide --potential-family option for calculations without POTCAR.', err=True)
            raise click.Abort()

    # Parse potential mapping if provided
    pot_mapping = None
    if potential_mapping:
        try:
            pot_mapping = json.loads(potential_mapping)
        except json.JSONDecodeError as e:
            click.echo(f'Error parsing potential mapping JSON: {e}', err=True)
            raise click.Abort()

    # Handle code - either load existing or create dummy
    vasp_code = None
    code_created = False

    if code:
        try:
            vasp_code = orm.load_code(code)
            if not quiet:
                click.echo(f'Using existing code: {vasp_code.label} (PK: {vasp_code.pk})')
        except NotExistent:
            click.echo(f'Error: Code "{code}" not found in database', err=True)
            raise click.Abort()
    else:
        # Create a dummy code for local imports only
        if is_remote:
            click.echo('Error: Code must be specified for remote imports', err=True)
            raise click.Abort()

        if not quiet:
            click.echo('No code specified, creating dummy code for local import...')
        try:
            vasp_code = orm.load_code('dummy-vasp-import@localhost')
        except NotExistent:
            vasp_code = _create_dummy_vasp_code()
            code_created = True

        if not quiet:
            click.echo(f'Created dummy code: {vasp_code.label} (PK: {vasp_code.pk})')

    # Prepare import parameters
    import_kwargs = {
        'code': vasp_code,
        'potential_family': potential_family,
        'potential_mapping': pot_mapping,
        'include_wavecar': include_wavecar,
        'include_chgcar': include_chgcar,
        'stdout_file_name': stdout_file,
    }

    import_kwargs['remote_path'] = str(calc_path.absolute())

    # Show summary and ask for confirmation
    if not yes and not quiet:
        click.echo('\n=== Import Summary ===')
        click.echo(f'Calculation path: {path}')
        click.echo(f'Code: {vasp_code.label} (PK: {vasp_code.pk})')
        if potential_family:
            click.echo(f'Potential family: {potential_family}')
        if pot_mapping:
            click.echo(f'Potential mapping: {pot_mapping}')
        if include_wavecar:
            click.echo('Will include WAVECAR if present')
        if include_chgcar:
            click.echo('Will include CHGCAR if present')
        if label:
            click.echo(f'Label: {label}')
        if description:
            click.echo(f'Description: {description}')
        if group:
            click.echo(f'Group: {group}')
        click.echo(f'Submit method: {"daemon" if submit_daemon else "direct"}')

        if not click.confirm('\nProceed with import?'):
            if code_created:
                # Clean up the dummy code if user cancels
                # Note: AiiDA nodes cannot be deleted once stored, so we just warn the user
                click.echo('Cancelled import. Note: dummy code was created and cannot be deleted.')
            else:
                click.echo('Cancelled import.')
            return

    try:
        # Perform the import
        if not quiet:
            click.echo('\nStarting import...')

        if submit_daemon:
            process = VaspCalcImporter.run_import_daemon(**import_kwargs)
            calc_node = process
            if not quiet:
                click.echo(f'Submitted import to daemon. Process PK: {process.pk}')
        else:
            calc_node = VaspCalcImporter.run_import(**import_kwargs)
            if not quiet:
                click.echo('Import completed successfully!')

        # Set metadata if provided
        if label:
            calc_node.label = label
        if description:
            calc_node.description = description

        # Add to group if specified
        if group:
            group_obj = orm.Group.collection.get_or_create(label=group)[0]
            group_obj.add_nodes([calc_node])
            if not quiet:
                click.echo(f'Added calculation to group: {group}')

        # Summarize what has been done
        if not quiet:
            click.echo('\n=== Import Results ===')
            click.echo(f'Calculation PK: {calc_node.pk}')
            click.echo(f'UUID: {calc_node.uuid}')
            click.echo(f'State: {calc_node.process_state}')

            if hasattr(calc_node, 'exit_status') and calc_node.exit_status is not None:
                click.echo(f'Exit status: {calc_node.exit_status}')

            # Show some basic info about inputs
            if 'structure' in calc_node.inputs:
                structure = calc_node.inputs.structure
                click.echo(f'Structure: {structure.get_formula()}')

            # Show outputs if available
            if calc_node.is_finished and 'misc' in calc_node.outputs:
                misc = calc_node.outputs.misc.get_dict()
                if 'total_energies' in misc:
                    energies = misc['total_energies']
                    if 'energy_extrapolated' in energies:
                        click.echo(f'Final energy: {energies["energy_extrapolated"]:.6f} eV')

            if code_created:
                click.echo(f'Note: Created dummy code {vasp_code.label} for this import')

    except Exception as e:
        click.echo(f'Error during import: {e}', err=True)
        raise click.Abort()


def _create_dummy_vasp_code():
    """Create a dummy VASP code for local imports."""
    from aiida import orm
    from aiida.common.exceptions import NotExistent

    # Try to get localhost computer
    try:
        computer = orm.Computer.collection.get(label='localhost')
    except NotExistent:
        # Create localhost computer if it doesn't exist
        computer = orm.Computer(
            label='localhost',
            hostname='localhost',
            transport_type='core.local',
            scheduler_type='core.direct',
            workdir='/tmp',
        )
        computer.store()
        computer.configure()

    # Create dummy code

    # Create a simple dummy script
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
        f.write('#!/bin/bash\necho "Dummy VASP code for import"\n')
        dummy_exec = f.name

    os.chmod(dummy_exec, 0o755)

    code = orm.InstalledCode(
        computer=computer,
        filepath_executable=dummy_exec,
    )
    code.label = 'dummy-vasp-import'
    code.description = 'Dummy VASP code created for importing calculations'
    code.default_calc_job_plugin = 'vasp.vasp'
    code.store()
    code.base.extras.set('is_dummy_import_code', True)

    return code
