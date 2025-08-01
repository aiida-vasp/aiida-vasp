"""
Provides aiida-vasp related tools as standalone commands.
"""

import json

import click

from . import cmd_aiida_vasp

# ruff: noqa: PLC0415


# @click.group(name='vasp-tools')  # pylint: disable=invalid-name
@cmd_aiida_vasp.group('launch')
def launch_cmd() -> None:
    """Top level command for aiida-vasp launch commands"""


def common_vasp_options(func):
    """Decorator to add common VASP calculation options."""
    func = click.option('--preset', '-p', default='VaspPreset', help='Preset to use for the calculation.')(func)
    func = click.option('--inputset', '-i', default='UCLRelaxSet', help='Input set to use for the calculation.')(func)
    func = click.option('--code', '-c', required=True, help='Code to use for the calculation.')(func)
    func = click.option(
        '--max-wallclock-seconds', '-m', type=int, default=None, help='Maximum wallclock time for the calculation.'
    )(func)
    func = click.option(
        '--num-machines', '-nm', type=int, default=None, help='Number of machines to use for the calculation.'
    )(func)
    func = click.option(
        '--tot-num-mpiprocs',
        '-np',
        type=int,
        default=None,
        help='Total number of MPI processes to use for the calculation.',
    )(func)
    func = click.option(
        '--resources', '-r', default=None, help='Resources for the calculation (JSON or key=value format).'
    )(func)
    func = click.option('--overrides', '-io', default=None, help='Path to a file containing input overrides')(func)
    func = click.option('--relax-settings', '-rs', default=None, help='Path to a file containing relaxation settings')(
        func
    )
    func = click.option('--band-settings', '-bs', default=None, help='Path to a file containing band settings')(func)
    func = click.option('--updates', '-u', default=None, help='Path to a file containing calls to set_xxx methods.')(
        func
    )
    func = click.option(
        '--structure', '-s', required=True, help='Path to a structure file to use for the calculation.'
    )(func)
    func = click.option('--group', '-g', default=None, help='Group to store the calculation in.')(func)
    func = click.option('--label', '-l', default=None, help='Label for the calculation.')(func)
    func = click.option('--description', '-d', default=None, help='Description for the calculation.')(func)
    func = click.option('--dryrun', is_flag=True, help='Show what would be done without actually submitting.')(func)
    func = click.option(
        '--run-directly', is_flag=True, help='Run the calculation directly in the current python process.'
    )(func)
    func = click.option(
        '--workchain-type',
        default='vasp',
        help='Type of workchain to launch.',
        type=click.Choice(['vasp', 'relax'], case_sensitive=False),
    )(func)
    return func


@cmd_aiida_vasp.command('launch-workchain')
@common_vasp_options
def launch_workchain(
    preset,
    inputset,
    code,
    max_wallclock_seconds,
    num_machines,
    resources,
    tot_num_mpiprocs,
    overrides,
    structure,
    group,
    label,
    description,
    dryrun,
    run_directly,
    workchain_type,
    relax_settings,
    band_settings,
    updates,
):
    """Launch a VASP calculation with the specified preset and input set."""
    from aiida_vasp.commands.utils import (
        apply_additional_updates,
        handle_calculation_submission,
        load_structure_from_file,
        process_dict_option,
        setup_calculation_options,
    )
    from aiida_vasp.common.builder_updater import (
        VaspBandUpdater,
        VaspBuilderUpdater,
        VaspConvUpdater,
        VaspHybridBandUpdater,
        VaspRelaxUpdater,
    )

    upd_cls_map = {
        'relax': VaspRelaxUpdater,
        'vasp': VaspBuilderUpdater,
        'band': VaspBandUpdater,
        'hybrid_band': VaspHybridBandUpdater,
        'conv': VaspConvUpdater,
    }
    try:
        # Load structure from file
        click.echo(f'Loading structure from: {structure}')
        structure_node = load_structure_from_file(structure)
        click.echo(f'Loaded structure: {structure_node.get_formula()}')

        # Initialize the builder updater
        click.echo(f'Initializing BuilderUpdater with preset: {preset}')
        upd_cls = upd_cls_map.get(workchain_type.lower(), VaspBuilderUpdater)
        upd = upd_cls(preset_name=preset, inputset_name=inputset)

        # Apply preset with structure
        upd.apply_preset(structure_node, code=code, label=label, overrides=process_dict_option(overrides))

        # Handle resource options
        options_dict = setup_calculation_options(resources, max_wallclock_seconds, num_machines, tot_num_mpiprocs)
        if options_dict:
            click.echo(f'Setting computational resources: {options_dict}')
            upd.set_options(**options_dict)

        # Apply any additional overrides
        apply_additional_updates(upd, process_dict_option(updates))

        if workchain_type.lower() == 'band':
            upd.set_band_settings(**process_dict_option(band_settings) if band_settings else {})
        if workchain_type.lower() == 'relax':
            upd.set_relax_settings(**process_dict_option(relax_settings) if relax_settings else {})

        # Set metadata
        if description:
            upd.builder.metadata.description = description

        if dryrun:
            click.echo(f'\n=== DRY RUN - Setup for {upd.builder._process_class.__name__} ===')
            click.echo(f'Code: {code}')
            click.echo(f'Structure: {structure_node.get_formula()} ({structure_node.label})')
            click.echo(f'Preset: {preset}')
            if inputset:
                click.echo(f'Input set: {inputset}')
            if label:
                click.echo(f'Label: {label}')
            if description:
                click.echo(f'Description: {description}')
            if options_dict:
                click.echo(f'Resources: {options_dict}')
            click.echo('Builder to be launched:')
            click.echo(pretty_print_builder(upd.builder))
            click.echo('=== END DRY RUN ===')
            return

        # Submit or run the calculation and handle groups
        handle_calculation_submission(upd, run_directly, group)

    except Exception as e:
        raise e
        click.echo(f'Error: {e}', err=True)
        raise click.Abort()


@cmd_aiida_vasp.command('list-presets')
@click.argument('preset', required=False, type=click.STRING)
def list_presets(preset):
    """List available presets for VASP calculations."""
    try:
        from aiida_vasp.common.builder_updater import list_presets

        preset_files = list_presets()
        if not preset_files:
            click.echo('No preset files found.')
            return

        if not preset:
            click.echo('\nAvailable presets:')
            click.echo('=' * 50)
        for preset_file in sorted(preset_files):
            name = preset_file.stem

            if preset:
                # Show a specific preset if provided
                if preset == name:
                    click.echo(preset_file.read_text())
                    return
                continue
            click.echo(f'• {name}: {preset_file}')
            click.echo('\nContent\n')
            click.echo('-' * 50)
            # Print the content of the file
            click.echo(preset_file.read_text())
            click.echo('-' * 50)

        click.echo('\nUse these preset names with the --preset option.')

    except Exception as e:
        click.echo(f'Error listing presets: {e}', err=True)
        raise click.Abort()


@cmd_aiida_vasp.command('list-inputsets')
@click.argument('inputset', required=False, type=click.STRING)
def list_inputset(inputset):
    """List available inputsets for VASP calculations."""
    try:
        from aiida_vasp.inputset.base import list_inputsets

        inputset_files = list_inputsets()
        if not inputset_files:
            click.echo('No inputset files found.')
            return

        if not inputset:
            click.echo('\nAvailable inputsets:')
            click.echo('=' * 50)
        for inputset_file in sorted(inputset_files):
            name = inputset_file.stem

            if inputset:
                # Show a specific inputset if provided
                if inputset == name:
                    click.echo(inputset_file.read_text())
                    return
                continue
            click.echo(f'• {name}: {inputset_file}')
            click.echo('\nContent\n')
            click.echo('-' * 50)
            # Print the content of the file
            click.echo(inputset_file.read_text())
            click.echo('-' * 50)

        click.echo('\nUse these inputsets names with the --inputset option.')

    except Exception as e:
        click.echo(f'Error listing inputset: {e}', err=True)
        raise click.Abort()


@cmd_aiida_vasp.command('status')
@click.argument('calculation_pk', type=int)
def status(calculation_pk):
    """Check the status of a VASP calculation."""
    try:
        from aiida.orm import load_node

        calc = load_node(calculation_pk)
        click.echo(f'Calculation PK: {calc.pk}')
        click.echo(f'UUID: {calc.uuid}')
        click.echo(f'Label: {calc.label}')
        click.echo(f'State: {calc.process_state}')

        if hasattr(calc, 'exit_status') and calc.exit_status is not None:
            click.echo(f'Exit status: {calc.exit_status}')

        if hasattr(calc, 'exit_message') and calc.exit_message:
            click.echo(f'Exit message: {calc.exit_message}')

        # Show creation and modification times
        click.echo(f'Created: {calc.ctime}')
        if hasattr(calc, 'mtime'):
            click.echo(f'Modified: {calc.mtime}')

        # Show inputs structure if available
        if 'structure' in calc.inputs:
            structure = calc.inputs.structure
            click.echo(f'Structure: {structure.get_formula()}')

        # Show some outputs if calculation is finished
        if calc.is_finished and 'misc' in calc.outputs:
            misc = calc.outputs.misc.get_dict()
            if 'total_energies' in misc:
                energies = misc['total_energies']
                if 'energy_extrapolated' in energies:
                    click.echo(f'Final energy: {energies["energy_extrapolated"]:.6f} eV')

    except Exception as e:
        click.echo(f'Error checking status: {e}', err=True)
        raise click.Abort()


def pretty_print_builder(builder) -> None:
    """
    Pretty print the builder object.

    Args:
        builder: The builder object to print.
        indent: Indentation level for pretty printing.
        stream: Output stream to write the pretty printed output.
    """
    import yaml
    from aiida.engine.processes.builder import PrettyEncoder

    return (
        f'Process class: {builder._process_class.__name__}\n'
        f'Inputs:\n{yaml.safe_dump(json.JSONDecoder().decode(PrettyEncoder().encode(builder)))}'
    )


class PrettyEncoder(json.JSONEncoder):
    """JSON encoder for returning a pretty representation of an AiiDA ``ProcessBuilder``."""

    def default(self, o):
        return dict(o)
