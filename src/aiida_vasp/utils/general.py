"""
General utils.

Contains general utils that is not directly coupled to the plugin or AiiDA.
"""

from __future__ import annotations

import os
import shutil
from typing import Callable


def copytree(
    src: str,
    dst: str,
    symlinks: bool = False,
    ignore: Callable[[str, list[str]], list[str]] | None = None,
) -> None:
    """
    Fixes annoying complaint about existing directory running tests

    In Python 3.8 there is a flag for shutil.copytree that handles this,
    but we still support 3.6 and 3.7.
    """
    for item in os.listdir(src):
        source = os.path.join(src, item)
        destination = os.path.join(dst, item)
        if os.path.isdir(source):
            shutil.copytree(source, destination, symlinks, ignore)
        else:
            shutil.copy2(source, destination)
