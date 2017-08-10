"""AiiDA Parser for aiida_vasp.Vasp2W90Calculation"""
from aiida.orm import DataFactory

from .vasp import VaspParser
from ..utils.io.win import WinParser


class Vasp2w90Parser(VaspParser):
    """Parse a finished aiida_vasp.Vasp2W90Calculation"""

    def parse_with_retrieved(self, retrieved):
        super(Vasp2w90Parser, self).parse_with_retrieved(retrieved)

        win_success, kpoints_node, param_node = self.parse_win()
        self.set_wannier_parameters(param_node)
        self.set_wannier_kpoints(kpoints_node)

        has_full_dat = self.has_full_dat()

        return self.result(success=win_success and has_full_dat)

    def parse_win(self):
        """Create the wannier90 .win file and kpoints output nodes."""
        win = self.get_file('wannier90.win')
        if not win:
            return False, None, None
        win_parser = WinParser(win)
        result = win_parser.result

        # remove kpoints block from parameters
        kpoints = result.pop('kpoints', None)
        success, kpoints_node = self.convert_kpoints(kpoints)

        # remove structure (cannot be given in parameters)
        result.pop('unit_cell_cart', None)
        result.pop('atoms_cart', None)

        param_node = DataFactory('parameter')(dict=result)
        return success, kpoints_node, param_node

    @staticmethod
    def convert_kpoints(kpoints):
        """Convert the k-points output from string to float."""
        if kpoints is None:
            return False, None
        kpoints_node = DataFactory('array.kpoints')()
        kpoints_node.set_kpoints([[float(x) for x in k.split()]
                                  for k in kpoints])
        return True, kpoints_node

    def set_wannier_parameters(self, node):
        """Add the Wannier parameters node."""
        if node:
            self.add_node('wannier_parameters', node)

    def set_wannier_kpoints(self, node):
        """Add the kpoints output node."""
        if node:
            self.add_node('wannier_kpoints', node)

    def has_full_dat(self):
        success = all(
            self.get_file('wannier90.' + ext) for ext in ['mmn', 'amn', 'eig'])
        return success
