"""
Representation of WAVECAR and WAVEDER objects.

--------------------------------------------
Wave function data node (stores WAVECAR and WAVEDER objects in the repository).
"""
# pylint: disable=abstract-method, too-few-public-methods, no-name-in-module
# explanation: pylint wrongly complains about (aiida) Node not implementing query
from aiida.orm import SinglefileData


class WavefunData(SinglefileData):
    pass
