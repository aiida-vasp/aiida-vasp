"""
Some convenience mixins
"""
from aiida.common.utils import classproperty

from .inputset.vaspsets import VASPInputSet


# pylint: disable=import-outside-toplevel, no-self-use
class WithVaspInputSet:
    """
    Mixins to attach a class property `vasp_inputset` which is the
    `VaspInputSet` class that is used for building the input files for VASP.
    """

    @classproperty
    def vasp_inputset(self):
        return VASPInputSet
