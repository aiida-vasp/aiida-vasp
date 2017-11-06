"""Provide pymatgen.io.vasp.outputs.Vasprun -> AiiDA nodes translation."""
from numpy import array
from pymatgen.electronic_structure.bandstructure import Spin

from aiida_vasp.utils.aiida_utils import dbenv


class VasprunToAiida(object):
    """Adapt from Vasprun object to AiiDA nodes."""

    def __init__(self, vasprun_obj, logger=None):
        self.vasprun_obj = vasprun_obj
        self.logger = logger

    @property
    def actual_kpoints(self):
        """List actually used kpoints as KpointsData."""
        kpoints = get_data_node('array.kpoints')
        kpoints.set_kpoints(
            self.vasprun_obj.actual_kpoints,
            weights=self.vasprun_obj.actual_kpoints_weights)
        return kpoints

    @property
    def last_structure(self):
        """Give the structure after or at the last recorded ionic step as StructureData."""
        last_structure = self.vasprun_obj.structures[-1]
        return self.final_structure or get_data_node(
            'structure', pymatgen=last_structure)

    @property
    def final_structure(self):
        """Give the final structure after all calculations as StructureData."""
        final_structure = getattr(self.vasprun_obj, 'final_structure', None)
        if not final_structure:
            return None
        return get_data_node('structure', pymatgen=final_structure)

    @property
    def forces(self):
        """Give forces in the last ioni step as ArrayData."""
        forces_array = get_data_node('array')
        forces_array.set_array('forces',
                               array(
                                   self.vasprun_obj.ionic_steps[-1]['forces']))
        return forces_array

    @property
    def output_dict(self):
        """Collect scalars and small arrays into a dictionary."""
        return {
            'stress': self.vasprun_obj.ionic_steps[-1]['stress'],
            'efermi': self.vasprun_obj.efermi,
            'energy': self.vasprun_obj.final_energy
        }

    @property
    def output_parameters(self):
        """Collect scalars and small arrays into a results ParameterData node."""
        output_params = get_data_node('parameter', dict=self.output_dict)
        return output_params

    @property
    def band_structure(self):
        """Give the band structure as BandsData node."""
        bands_node = get_data_node('array.bands')
        bands_node.set_kpointsdata(self.actual_kpoints)  # has to be set first
        try:
            bands_object = self.vasprun_obj.get_band_structure()
            structure = get_data_node(
                'structure', pymatgen=bands_object.structure)
            bands_node.set_cell_from_structure(structure)
            bands_data = bands_object.bands
            bands_node_data = []
            for spin in [Spin.up, Spin.down]:
                if spin in bands_data:
                    bands_node_data.append(bands_data[spin].transpose())
            bands_node_data = array(bands_node_data)
            bands_node.set_bands(bands=bands_node_data)
        except AttributeError:
            if self.logger:
                self.logger.warning(
                    'Band structure could not be parsed, possibly because the final structure was missing from the xml'
                )
        return bands_node

    @property
    def tdos(self):
        """Give the total density of states as an ArrayNode."""
        tdos = self.vasprun_obj.tdos
        tdos_node = get_data_node('array.xy')
        tdos_node.set_x(tdos.energies, 'dos_energy', 'eV')
        spins = [(Spin.up, 'dos_spin_up'), (Spin.down, 'dos_spin_down')]
        y_arrays = []
        y_names = []
        y_units = []
        for spin, y_label in spins:
            if spin in tdos.densities:
                y_arrays.append(tdos.densities[spin])
                y_names.append(y_label)
                y_units.append('states/eV')
        tdos_node.set_y(y_arrays, y_names, y_units)
        return tdos_node


@dbenv
def get_data_node(data_type, *args, **kwargs):
    from aiida.orm import DataFactory
    return DataFactory(data_type)(*args, **kwargs)
