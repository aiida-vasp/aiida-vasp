""" # noqa: D205
Representation of WAVECAR and WAVEDER files
-------------------------------------------
Wave function data node (stores WAVECAR and WAVEDER files).
"""
# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
from aiida.orm import SinglefileData


class WavefunData(SinglefileData):
    pass
