from aiida.orm import DataFactory

from .vasp import VaspParser
from ..utils.io.win import WinParser


class Vasp2w90Parser(VaspParser):
    def parse_with_retrieved(self, retrieved):
        super(Vasp2w90Parser, self).parse_with_retrieved(retrieved)

        has_win, win_node = self.get_win_node()
        has_full_dat, dat_node = self.get_wdat_node()
        self.set_win(win_node)
        self.set_wdat(dat_node)

        return self.result(success=has_win and has_full_dat)

    def get_win_node(self):
        win = self.get_file('wannier90.win')
        if not win:
            return False, None
        wp = WinParser(win)
        winnode = DataFactory('parameter')(dict=wp.result)
        return True, winnode

    def get_wdat_node(self):
        wdatnode = DataFactory('vasp.archive')()
        success = True
        for ext in ['mmn', 'amn', 'eig']:
            wfile = self.get_file('wannier90.' + ext)
            if wfile:
                wdatnode.add_file(wfile)
            else:
                success = False
        return success, wdatnode

    def set_win(self, node):
        if node:
            self.add_node('wannier_parameters', node)

    def set_wdat(self, node):
        if node:
            self.add_node('wannier_data', node)
