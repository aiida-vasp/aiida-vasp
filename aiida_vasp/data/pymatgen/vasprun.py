"""pymatgen.io.vasp.outputs.Vasprun -> AiiDA nodes translation"""
from numpy import array

from aiida_vasp.utils.aiida_utils import dbenv


class VasprunToAiida(object):
    """Adapter from Vasprun object to AiiDA nodes"""

    def __init__(self, vasprun_obj):
        self.vasprun_obj = vasprun_obj

    @property
    def actual_kpoints(self):
        """
        List of actually used kpoints as KpointsData
        """
        kpoints = get_data_node('array.kpoints')
        kpoints.set_kpoints(
            self.vasprun_obj.actual_kpoints,
            weights=self.vasprun_obj.actual_kpoints_weights)
        return kpoints

    @property
    def last_structure(self):
        """
        The structure after or at the last recorded ionic step as StructureData
        """
        last_structure = self.vasprun_obj.structures[-1]
        return self.final_structure or get_data_node(
            'structure', pymatgen=last_structure)

    @property
    def final_structure(self):
        """
        The final structure after all calculations as StructureData
        """
        final_structure = getattr(self.vasprun_obj, 'final_structure', None)
        if not final_structure:
            return None
        return get_data_node('structure', pymatgen=final_structure)

    @property
    def forces(self):
        """Forces in the last ioni step as ArrayData"""
        forces_array = get_data_node('array')
        forces_array.set_array('forces',
                               array(
                                   self.vasprun_obj.ionic_steps[-1]['forces']))
        return forces_array


@dbenv
def get_data_node(data_type, *args, **kwargs):
    from aiida.orm import DataFactory
    return DataFactory(data_type)(*args, **kwargs)
