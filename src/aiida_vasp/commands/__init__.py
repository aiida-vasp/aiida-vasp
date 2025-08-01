"""
Root command for aiida-vasp
"""

# ruff: noqa: E402
import click
from aiida.cmdline.groups import VerdiCommandGroup
from aiida.cmdline.params import options, types


@click.group(
    'aiida-vasp',
    cls=VerdiCommandGroup,
    help='AiiDA VASP command line tools',
    context_settings={'help_option_names': ['-h', '--help']},
)
@options.PROFILE(type=types.ProfileParamType(load_profile=True), expose_value=False)
def cmd_aiida_vasp() -> None:
    """
    AiiDA VASP command line tools.

    This is the root command for all aiida-vasp related commands.
    """
    pass


# These lines simply triggers python to include the sub-commands
from .launch import launch_workchain as launch_workchain
