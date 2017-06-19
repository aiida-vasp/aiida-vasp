from aiida.orm import DataFactory

from .wannier import WannierParser

class WswannierParser(WannierParser):
    def get_hr_node(self):
        hr_filename = 'wannier90_hr.dat'
        hr_dat = self.get_file(hr_filename)
        if not hr_dat:
            return None
        tbnode = DataFactory('folder')()
        tbnode.add_path(hr_dat, hr_filename)
        # adding optional files
        for filename in [
                'wannier90_centres.xyz',
                'wannier90_wsvec.dat',
                'wannier90_tb.dat',
                'wannier90.win'
        ]:
            f_dat = self.get_file(filename)
            if f_dat:
                tbnode.add_path(f_dat, filename)
        return tbnode
