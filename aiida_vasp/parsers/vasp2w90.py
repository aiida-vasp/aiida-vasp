from aiida_vasp.parsers.vasp5 import Vasp5Parser
from aiida_vasp.utils.io.win import WinParser


class Vasp2w90Parser(Vasp5Parser):
    def parse_with_retrieved(self, retrieved):
        super(Vasp2w90Parser, self).parse_with_retrieved(retrieved)

        self.set_win(self.get_win_node())
        self.set_wdat(self.get_wdat_node())

        return self.result(success=True)

    def get_win_node(self):
        if self._calc.get_inputs_dict().get('wannier_settings'):
            return None
        win = self.get_file('wannier90.win')
        if not win:
            return None
        wp = WinParser(win)
        winnode = self._calc.new_wannier_settings(dict=wp.result)
        return winnode

    def get_wdat_node(self):
        if self._calc.get_inputs_dict().get('wannier_settings'):
            return None
        wdatnode = self._calc.new_wannier_data()
        for ext in ['mmn', 'amn', 'eig']:
            wfile = self.get_file('wannier90.' + ext)
            if wfile:
                wdatnode.add_file(wfile)
        return wdatnode

    def set_win(self, node):
        if node:
            self.add_node('wannier_settings', node)

    def set_wdat(self, node):
        if node:
            self.add_node('wannier_data', node)
