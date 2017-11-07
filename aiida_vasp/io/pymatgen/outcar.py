"""Translate from pymatgen Outcar class to Aiida datastructures."""
from aiida_vasp.io.pymatgen.vasprun import get_data_node


class OutcarToAiida(object):
    """Adapt from Outcar object to AiiDA nodes."""

    def __init__(self, outcar_obj):
        self.outcar_obj = outcar_obj

    @property
    def output_dict(self):
        """Collect scalars and small arrays into a dictionary."""
        output_dict = {}
        if self.outcar_obj.lepsilon:
            output_dict[
                'dielectric tensor'] = self.outcar_obj.dielectric_tensor
        return output_dict

    @property
    def born_charges(self):
        """Give the born charges as an ArrayData node."""
        born_node = get_data_node('array')
        born_node.set_array('born_charges', self.outcar_obj.born)
        return born_node
