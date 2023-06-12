"""
Representation of CHGCAR objects.

-------------------------------
Charge density data node (stores CHGCAR objects in the repository).
"""
# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
from aiida.plugins import DataFactory

SinglefileData = DataFactory('core.singlefile')


class ChargedensityData(SinglefileData):  # pylint: disable=too-few-public-methods
    pass
