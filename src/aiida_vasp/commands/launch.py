"""
Provides aiida-vasp related tools as standalone commands.
"""

import json

import click

from . import cmd_aiida_vasp

# ruff: noqa: PLC0415


def common_vasp_options(func):
    """Decorator to add common VASP calculation options."""
    func = click.option('--preset', '-p', default='default', help='Preset to use for the calculation.')(func)
    func = click.option('--protocol', '-i', default='balanced', help='Input set to use for the calculation.')(func)
    func = click.option('--code', '-c', required=False, help='Code to use for the calculation.')(func)
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
        '--options', '-op', default=None, help='Options for the calculation (JSON or key=value format).'
    )(func)
    func = click.option(
        '--resources', '-r', default=None, help='Options for the calculation (JSON or key=value format).'
    )(func)
    func = click.option('--overrides', '-io', default=None, help='Path to a file containing input overrides')(func)
    func = click.option('--relax-settings', '-rs', default=None, help='Path to a file containing relaxation settings')(
        func
    )
    func = click.option(
        '--incar-overrides', help='Additional incar overrides to be passed as set_incar method of the InputGenerator.'
    )
    func = click.option('--band-settings', '-bs', default=None, help='Path to a file containing band settings')(func)
    func = click.option('--updates', '-u', default=None, help='Path to a file containing calls to set_xxx methods.')(
        func
    )
    func = click.option(
        '--structure', '-s', required=False, help='Path to a structure file to use for the calculation.'
    )(func)
    func = click.option(
        '--from-vasp-folder', '-fvf', required=False, help='Path to existing VASP folder to use as input template.'
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
    protocol,
    code,
    max_wallclock_seconds,
    num_machines,
    resources,
    options,
    tot_num_mpiprocs,
    overrides,
    structure,
    from_vasp_folder,
    incar_overrides,
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
    from aiida_vasp.protocols.generator import (
        VaspBandUpdater,
        VaspConvergenceUpdater,
        VaspHybridBandsUpdater,
        VaspRelaxUpdater,
        VaspUpdater,
    )
    from aiida_vasp.utils.dict_merge import recursive_merge

    upd_cls_map = {
        'relax': VaspRelaxUpdater,
        'vasp': VaspUpdater,
        'band': VaspBandUpdater,
        'hybrid_band': VaspHybridBandsUpdater,
        'conv': VaspConvergenceUpdater,
    }
    try:
        # Validate input sources
        if not structure and not from_vasp_folder:
            click.echo('Error: Either --structure or --from-vasp-folder must be specified', err=True)
            raise click.Abort()

        if structure and from_vasp_folder:
            click.echo('Error: Cannot specify both --structure and --from-vasp-folder', err=True)
            raise click.Abort()

        # Load structure from file or VASP folder
        if from_vasp_folder:
            structure_node, overrides_map = load_inputs_from_vasp_folder(from_vasp_folder)
            click.echo(f'Loaded structure from VASP folder: {from_vasp_folder}')
            click.echo(f'Structure: {structure_node.get_formula()}')
            wc_type = workchain_type.lower()
            if wc_type in overrides_map:
                vasp_overrides = overrides_map[wc_type]
            else:
                vasp_overrides = {}
        else:
            structure_node = load_structure_from_file(structure)
            click.echo(f'Loading structure from: {structure}')
            click.echo(f'Loaded structure: {structure_node.get_formula()}')
            vasp_overrides = {}

        # Initialize the builder updater
        click.echo(f'Initializing BuilderUpdater with preset: {preset}')
        upd_cls = upd_cls_map.get(workchain_type.lower(), VaspUpdater)
        upd = upd_cls(preset_name=preset, protocol=protocol)

        # Merge overrides from VASP folder with user-specified overrides
        final_overrides = recursive_merge(vasp_overrides, process_dict_option(overrides))

        # Apply preset with structure
        upd.get_builder(structure=structure_node, code=code, overrides=final_overrides)
        upd.set_label(label)

        # Handle resource options
        options_dict = setup_calculation_options(
            options, resources, max_wallclock_seconds, num_machines, tot_num_mpiprocs
        )
        if options_dict:
            click.echo(f'Setting computational resources: {options_dict}')
            upd.set_options(options_dict)

        # Apply any additional overrides
        apply_additional_updates(upd, process_dict_option(updates))

        if workchain_type.lower() == 'band':
            upd.set_band_settings(process_dict_option(band_settings) if band_settings else {})
        if workchain_type.lower() == 'relax':
            upd.set_relax_settings(process_dict_option(relax_settings) if relax_settings else {})

        # Set metadata
        if description:
            upd.builder.metadata.description = description

        # Apply incar overrides
        if incar_overrides is not None:
            upd.set_incar(process_dict_option(incar_overrides))

        if dryrun:
            click.echo(f'\n=== DRY RUN - Setup for {upd.builder._process_class.__name__} ===')
            click.echo(f'Code: {code}')
            click.echo(f'Structure: {structure_node.get_formula()} ({structure_node.label})')
            if from_vasp_folder:
                click.echo(f'VASP folder: {from_vasp_folder}')
                incar_params = vasp_overrides.get('base', {}).get('vasp', {}).get('parameters', {}).get('incar', {})
                click.echo(f'INCAR parameters loaded: {len(incar_params)}')
            click.echo(f'Preset: {preset}')
            if protocol:
                click.echo(f'Protocol: {protocol}')
            if label:
                click.echo(f'Label: {label}')
            if description:
                click.echo(f'Description: {description}')
            if final_overrides:
                click.echo(f'Overrides applied: {len(final_overrides)} top-level keys')
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
@click.option('--show-content', default=False, is_flag=True, help='Include the content of the protocol files.')
def list_presets(preset, show_content):
    """List available presets for VASP calculations."""
    from aiida_vasp.protocols.generator import list_protocol_presets

    preset_files = list_protocol_presets()
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
        if show_content:
            click.echo('\nContent\n')
            click.echo('-' * 50)
            # Print the content of the file
            click.echo(preset_file.read_text())
            click.echo('-' * 50)

    click.echo('\nHint: Use these preset names with the --preset option.')


@cmd_aiida_vasp.command('list-protocols')
@click.argument('workflow-tag', required=False, type=click.STRING)
@click.option('--show-content', default=False, is_flag=True, help='Include the content of the protocol files.')
def list_protocols(workflow_tag, show_content):
    """List available protocols for VASP calculations."""
    from yaml import safe_load

    from aiida_vasp.protocols import ProtocolMixin

    protocol_files = ProtocolMixin.list_protocol_files(protocol_tag=workflow_tag)

    if protocol_files:
        click.echo('\nAvailable files containing protocols:')
        click.echo('=' * 80)
    else:
        click.echo(f'No protocol files found for {workflow_tag}')
    for _alias, tag, path in protocol_files:
        alias = _alias or 'default'
        click.echo(f'• workflow {tag:5s} -> protocol alias {alias:10s}: {path}')
        with open(path, 'r') as f:
            click.echo(f'  - available protocols: {list(safe_load(f).get("protocols"))} ')
        if show_content:
            click.echo('=' * 80)
            click.echo(f'\nContent of {path}\n')
            click.echo('-' * 80)
            # Print the content of the file
            click.echo(path.read_text())
            click.echo('-' * 80)

    click.echo('\nHint: Use these protocol names with the --protocol option for launching calculations.')


# TODO - print a tree-like diagram for calculation
@cmd_aiida_vasp.command('status')
@click.argument('process_pk')
def status(process_pk):
    """Check the status of a VaspCalculation or VasoWorkChain."""
    from aiida import orm

    def print_calculation_info(calculation_pk):
        calc = orm.load_node(calculation_pk)
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

    node = orm.load_node(process_pk)
    if isinstance(node, orm.CalcJobNode):
        print_calculation_info(node.pk)
    else:
        for node in node.called_descendants:
            if isinstance(node, orm.CalcJobNode):
                print_calculation_info(node.pk)
                click.echo('-' * 40)


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


def load_inputs_from_vasp_folder(folder_path):
    """
    Load structure and parameters from an existing VASP calculation folder.

    :param folder_path: Path to the VASP calculation folder
    :return: Tuple of (structure_node, overrides_dict)
    """
    from pathlib import Path

    from aiida_vasp.calcs.immigrant import get_incar_input, get_kpoints_input, get_poscar_input, get_potcar_input

    folder = Path(folder_path)

    if not folder.exists():
        raise click.ClickException(f'VASP folder not found: {folder_path}')

    # Check for required files
    required_files = ['INCAR', 'POSCAR']
    missing_files = [f for f in required_files if not (folder / f).exists()]

    if missing_files:
        raise click.ClickException(f'Missing required files in VASP folder: {missing_files}')

    # Load structure from POSCAR
    try:
        structure_node = get_poscar_input(folder)
        click.echo(f'Loaded structure with formula: {structure_node.get_formula()}')
    except Exception as e:
        raise click.ClickException(f'Error reading POSCAR: {e}')

    # Load INCAR parameters
    try:
        incar_node = get_incar_input(folder)
        incar_dict = incar_node.get_dict()
        click.echo(f'Loaded INCAR with {len(incar_dict)} parameters')
    except Exception as e:
        raise click.ClickException(f'Error reading INCAR: {e}')

    # Prepare overrides for each workchain type
    overrides_map = {}

    # KPOINTS
    kpoints_file = folder / 'KPOINTS'
    kpoints_node = None
    if kpoints_file.exists():
        try:
            kpoints_node = get_kpoints_input(folder, structure_node)
        except Exception as e:
            click.echo(f'Warning: Could not process KPOINTS file: {e}')

    # POTCARS
    if (folder / 'POTCAR').exists():
        potcars = get_potcar_input(folder, structure_node)
    else:
        potcars = {}

    # Compose overrides for each workchain type
    # VaspWorkChain
    vasp_override = {'parameters': {'incar': incar_dict}, 'kpoints': kpoints_node, 'potential': potcars}
    overrides_map['vasp'] = vasp_override

    # VaspRelaxWorkChain
    relax_override = {'vasp': {'parameters': {'incar': incar_dict}, 'kpoints': kpoints_node, 'potential': potcars}}
    overrides_map['relax'] = relax_override

    # VaspBandsWorkChain
    band_override = {'scf': {'parameters': {'incar': incar_dict}, 'kpoints': kpoints_node, 'potential': potcars}}
    overrides_map['band'] = band_override

    # VaspConvergeWorkChain
    conv_override = {'parameters': {'incar': incar_dict, 'kpoints': kpoints_node, 'potential': potcars}}
    overrides_map['conv'] = conv_override

    return structure_node, overrides_map
