"""
Root command for aiida-vasp
"""

# ruff: noqa: E402
import click
from aiida.cmdline.groups import VerdiCommandGroup
from aiida.cmdline.params.options import PROFILE
from aiida.cmdline.params.types.profile import ProfileParamType


@click.group(
    'aiida-vasp',
    cls=VerdiCommandGroup,
    help='AiiDA VASP command line tools',
    context_settings={'help_option_names': ['-h', '--help']},
)
@PROFILE(type=ProfileParamType(load_profile=True), expose_value=False)
def cmd_aiida_vasp() -> None:
    """
    AiiDA VASP command line tools.

    This is the root command for all aiida-vasp related commands.
    """
    pass


# These lines simply triggers python to include the sub-commands
from .immigrant import import_calc
from .launch import launch_workchain
from .potcar import potcar
from .tools import export

__all__ = ['cmd_aiida_vasp', 'export', 'import_calc', 'launch_workchain', 'potcar']
