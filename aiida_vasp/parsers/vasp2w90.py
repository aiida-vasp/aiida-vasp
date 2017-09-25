"""AiiDA Parser for aiida_vasp.Vasp2W90Calculation"""
from aiida.orm import DataFactory
from aiida.orm.data.base import List

from .vasp import VaspParser
from ..utils.io.win import WinParser


class Vasp2w90Parser(VaspParser):
    """Parse a finished aiida_vasp.Vasp2W90Calculation"""

    def parse_with_retrieved(self, retrieved):
        super(Vasp2w90Parser, self).parse_with_retrieved(retrieved)

        win_success, kpoints_node, param_node, proj_node = self.parse_win()
        self.set_node('wannier_parameters', param_node)
        self.set_node('wannier_kpoints', kpoints_node)
        self.set_node('wannier_projections', proj_node)

        has_full_dat = self.has_full_dat()

        return self.result(success=win_success and has_full_dat)

    def parse_win(self):
        """Create the wannier90 .win file and kpoints output nodes."""
        win = self.get_file('wannier90.win')
        if not win:
            return False, None, None, None
        win_result = WinParser(win).result

        # remove kpoints block from parameters
        kpoints = win_result.pop('kpoints', None)
        success, kpoints_node = self.convert_kpoints(kpoints)

        # remove structure (cannot be given in parameters)
        win_result.pop('unit_cell_cart', None)
        win_result.pop('atoms_cart', None)

        projections = win_result.pop('projections', None)
        if projections is None:
            proj_node = None
        else:
            proj_node = List()
            proj_node.extend(projections)

        param_node = DataFactory('parameter')(dict=win_result)
        return success, kpoints_node, param_node, proj_node

    @staticmethod
    def convert_kpoints(kpoints):
        """Convert the k-points output from string to float."""
        if kpoints is None:
            return False, None
        kpoints_node = DataFactory('array.kpoints')()
        kpoints_node.set_kpoints([[float(x) for x in k.split()]
                                  for k in kpoints])
        return True, kpoints_node

    def set_node(self, name, node):
        """Add a node if it is not None"""
        if node is not None:
            self.add_node(name, node)

    def has_full_dat(self):
        success = all(
            self.get_file('wannier90.' + ext) for ext in ['mmn', 'amn', 'eig'])
        return success
