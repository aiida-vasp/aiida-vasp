"""AiiDA Parser for aiida_vasp.WannierCalculation"""
import re

import numpy as np
from aiida.orm import DataFactory

from aiida_vasp.parsers.base import BaseParser
from aiida_vasp.utils.io.win import WinParser
from aiida_vasp.utils.io import parser


class WannierBase(BaseParser):
    """Parse a finished aiida_vasp.WannierCalculation"""

    def parse_with_retrieved(self, retrieved):
        super(WannierBase, self).parse_with_retrieved(retrieved)

        return self.result(success=True)

    def get_win_node(self):
        """Create the output settings node"""
        if self._calc.get_inputs_dict().get('wannier_settings'):
            return None
        win = self.get_file('wannier90.win')
        if not win:
            return None
        win_parser = WinParser(win)
        winnode = self._calc.new_wannier_settings(dict=win_parser.result)
        return winnode

    def get_wdat_node(self):
        """Create the wannier data output node comprised of .mmn, .amn, .eig files"""
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


class WannierParser(WannierBase, parser.BaseParser):
    """Parse a finished aiida_vasp.WannierCalculation"""

    def parse_with_retrieved(self, retrieved):
        super(WannierParser, self).parse_with_retrieved(retrieved)
        bands_node = self.get_bands_node()
        if bands_node:
            self.add_node('bands', bands_node)
        self.add_node('tb_model', self.get_hr_node())
        return self.result(success=True)

    def get_hr_node(self):
        """Create the hamiltonian output node"""
        dat = self.get_file('wannier90_hr.dat')
        if not dat:
            return None
        hnode = DataFactory('singlefile')()
        hnode.add_path(dat)
        return hnode

    def get_bands_plot(self):
        bands_file = self.get_file('wannier90_band.dat')
        bands_kp_file = self.get_file('wannier90_band.kpt')
        bands_gnu_file = self.get_file('wannier90_band.gnu')
        return bands_file, bands_kp_file, bands_gnu_file

    def get_bands_node(self):
        """Create the bandstructure output node"""
        bands_node = DataFactory('array.bands')()
        bands_file, bands_kp_file, _ = self.get_bands_plot()
        if not (bands_file and bands_kp_file):
            return None
        with open(bands_kp_file) as bands_kp_f:
            self.line(bands_kp_f)  # contains num kpoints
            kpoints = re.split(self.empty_line, bands_kp_f.read())
            kpoints = filter(None, kpoints)
            kpoints = map(self.splitlines, kpoints)
            kpoints = filter(None, kpoints[0])
            kpoints = np.array(kpoints)
            bands_node.set_kpoints(kpoints[:, :3], weights=kpoints[:, 3])
        with open(bands_file) as bands_f:
            data = re.split(self.empty_line, bands_f.read())
            data = filter(None, data)
            data = map(self.splitlines, data)
            data = np.array(data)
        bands_node.set_bands(data[:, :, 1].transpose())
        kppath = self._calc.inp.settings.get_dict().get('kpoint_path')
        kpl = [[kpp[0], kpp[1:4]] for kpp in kppath]
        kpl.append([kppath[-1][4], kppath[-1][5:8]])
        counter = {i[0]: 0 for i in kpl}
        kplab = []
        for kpi in kpl:
            ctr_i = counter[kpi[0]]
            idx = self._find_special_kpoint(kpoints[:, :3], kpi[1], num=ctr_i)
            kplab.append((idx, kpi[0]))
            counter[kpi[0]] += 1
        bands_node.labels = kplab
        return bands_node

    @staticmethod
    def _find_special_kpoint(kpoints, special_kp, num=0):
        """Find the special kpoints in the kpoints list"""
        res = []
        i_x, i_y = np.where(kpoints == special_kp)
        for i in range(len(i_x) - 2):
            if i_x[i] == i_x[i + 1] == i_x[i + 2]:
                if np.all(i_y[i:i + 3] == [0, 1, 2]):
                    res.append(i_x[i])
        return res[num]
