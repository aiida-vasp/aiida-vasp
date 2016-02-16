from aiida.parsers.plugins.vasp.vasp5 import Vasp5Parser
from aiida.tools.codespecific.vasp.io.win import WinParser


class Vasp2W90Parser(Vasp5Parser):
    def parse_with_retrieved(self, retrieved):
        super(Vasp2W90Parser, self).parse_with_retrieved(retrieved)

        winnode = self.get_win_node()
        self.set_win(winnode)

        return self.result(success=True)

    def get_win_node(self):
        if self._calc.get_inputs_dict().get('wannier_settings'):
            return None
        win = self.get_file('wannier90.win')
        wp = WinParser(win)
        winnode = self._calc.new_wannier_settings(dict=wp.result)
        return winnode

    def set_win(self, node):
        if node:
            self.add_node('wannier_settings', node)
