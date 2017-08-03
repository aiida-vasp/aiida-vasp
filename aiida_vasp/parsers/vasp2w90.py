"""AiiDA Parser for aiida_vasp.Vasp2W90Calculation"""
from aiida.orm import DataFactory

from .vasp import VaspParser
from ..utils.io.win import WinParser


class Vasp2w90Parser(VaspParser):
    """Parse a finished aiida_vasp.Vasp2W90Calculation"""

    def parse_with_retrieved(self, retrieved):
        super(Vasp2w90Parser, self).parse_with_retrieved(retrieved)

        has_win, win_node = self.get_win_node()
        has_full_dat = self.has_full_dat()
        self.set_win(win_node)

        return self.result(success=has_win and has_full_dat)

    def get_win_node(self):
        """Create the wannier90 .win file output node"""
        win = self.get_file('wannier90.win')
        if not win:
            return False, None
        win_parser = WinParser(win)
        winnode = DataFactory('parameter')(dict=win_parser.result)
        return True, winnode

    def has_full_dat(self):
        success = all(
            self.get_file('wannier90.' + ext) for ext in ['mmn', 'amn', 'eig'])
        return success

    def set_win(self, node):
        if node:
            self.add_node('wannier_parameters', node)
