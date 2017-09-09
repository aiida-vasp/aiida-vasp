"""AiiDA Parser for aiida_vasp.Vasp5Calculation"""
from aiida_vasp.parsers.vasp5 import Vasp5Parser
from aiida_vasp.utils.io.win import WinParser


class Vasp2w90Parser(Vasp5Parser):
    """Parse a finished aiida_vasp.Vasp5Calculation"""

    def parse_with_retrieved(self, retrieved):
        super(Vasp2w90Parser, self).parse_with_retrieved(retrieved)

        self.set_win(self.get_win_node())
        self.set_wdat(self.get_wdat_node())

        return self.result(success=True)

    def get_win_node(self):
        """Create the wannier90 .win file output node"""
        if self._calc.get_inputs_dict().get('wannier_parameters'):
            return None
        win = self.get_file('wannier90.win')
        if not win:
            return None
        win_parser = WinParser(win)
        winnode = self._calc.new_wannier_parameters(dict=win_parser.result)
        return winnode

    def get_wdat_node(self):
        """Create the wannier90 data output node comprised of the .mmn, .amn, .eig files"""
        if self._calc.get_inputs_dict().get('wannier_parameters'):
            return None
        wdatnode = self._calc.new_wannier_data()
        for ext in ['mmn', 'amn', 'eig']:
            wfile = self.get_file('wannier90.' + ext)
            if wfile:
                wdatnode.add_file(wfile)
        return wdatnode

    def set_win(self, node):
        if node:
            self.add_node('wannier_parameters', node)

    def set_wdat(self, node):
        if node:
            self.add_node('wannier_data', node)
