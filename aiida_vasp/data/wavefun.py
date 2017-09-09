"""Charge density data node (stores WAVECAR files)"""
# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
from aiida.orm.data.singlefile import SinglefileData


class WavefunData(SinglefileData):
    pass
