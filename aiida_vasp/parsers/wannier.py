from aiida_vasp.parsers.base import BaseParser
from aiida_vasp.utils.io.win import WinParser
from aiida_vasp.utils.io import parser
from aiida.orm import DataFactory
import re
import numpy as np


class WannierBase(BaseParser):
    def parse_with_retrieved(self, retrieved):
        super(WannierBase, self).parse_with_retrieved(retrieved)

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
            wfile = self.get_file('wannier90.'+ext)
            if wfile:
                wdatnode.add_file(wfile)
        return wdatnode

    def set_win(self, node):
        if node:
            self.add_node('wannier_settings', node)

    def set_wdat(self, node):
        if node:
            self.add_node('wannier_data', node)


class WannierParser(WannierBase, parser.BaseParser):
    def parse_with_retrieved(self, retrieved):
        super(WannierBase, self).parse_with_retrieved(retrieved)
        self.add_node('bands', self.get_bands_node())
        self.add_node('tb_model', self.get_hr_node())
        return self.result(success=True)

    def get_hr_node(self):
        dat = self.get_file('wannier90_hr.dat')
        if not dat:
            return None
        hnode = DataFactory('singlefile')()
        hnode.add_path(dat)
        return hnode

    def get_bands_plot(self):
        b = self.get_file('wannier90_band.dat')
        bkp = self.get_file('wannier90_band.kpt')
        bgnu = self.get_file('wannier90_band.gnu')
        return b, bkp, bgnu

    def get_bands_node(self):
        bnode = DataFactory('array.bands')()
        bdat, bkp, bgnu = self.get_bands_plot()
        if not (bdat and bkp):
            return None
        with open(bkp) as bk:
            self.line(bk)  # contains num kpoints
            kp = re.split(self.empty_line, bk.read())
            kp = filter(None, kp)
            kp = map(self.splitlines, kp)
            kp = filter(None, kp[0])
            kp = np.array(kp)
            bnode.set_kpoints(kp[:, :3], weights=kp[:, 3])
        with open(bdat) as bd:
            data = re.split(self.empty_line, bd.read())
            data = filter(None, data)
            data = map(self.splitlines, data)
            data = np.array(data)
        bnode.set_bands(data[:, :, 1].transpose())
        kppath = self._calc.inp.settings.get_dict().get('kpoint_path')
        kpl = [[kpp[0], kpp[1:4]] for kpp in kppath]
        kpl.append([kppath[-1][4], kppath[-1][5:8]])
        counter = {i[0]: 0 for i in kpl}
        kplab = []
        for kpi in kpl:
            ci = counter[kpi[0]]
            idx = self._find_special_kpoint(kp[:, :3], kpi[1], num=ci)
            kplab.append((idx, kpi[0]))
            counter[kpi[0]] += 1
        bnode.labels = kplab
        return bnode

    def _find_special_kpoint(self, kp, sp, num=0):
        res = []
        ix, iy = np.where(kp == sp)
        for i in range(len(ix)-2):
            if ix[i] == ix[i + 1] == ix[i + 2]:
                if np.all(iy[i:i+3] == [0, 1, 2]):
                    res.append(ix[i])
        return res[num]
