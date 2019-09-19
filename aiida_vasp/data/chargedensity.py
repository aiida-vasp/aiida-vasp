""" # noqa: D205
Representation of CHGCAR files
------------------------------

Charge density data node (stores CHGCAR files).
"""
# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
from aiida.orm import SinglefileData


class ChargedensityData(SinglefileData):
    pass
