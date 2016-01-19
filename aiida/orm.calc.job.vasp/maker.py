froma aiida.orm import CalculationFactory

class VaspMaker(object):
    def __init__(self, calcname, *args, **kwargs):
        self.calc_cls = CalculationFactory(calcname)
        self.label = kwargs.get('label')
        self.cif = kwargs.get('cif')
        self.settings = kwargs.get('incar', {})
        self.paws = kwargs.get('paws', {})
        self.struct = kwargs.get('structure', {})
        self.computer = kwargs.get('computer')
        self.code = kwargs.get('code')
        self.make_structure()

    def make_structure(self):
        if self.cif:
            from aiida.tools.codespecific.vasp.io import cif
            from os.path import basename
            cifnode = cif.cif_from_file(self.cif)
            nmcifs = cif.get_cifs_with_name(basename(self.cif))
            eqcifs = cif.filter_ifs_for_structure(nmcifs, cifnode.get_ase())
            if eqcifs:
                cifnode = eqcifs[0]
            self.struct = cif.cif_to_structure(cifnode=cifnode)

    def verify_incar(self):
        if not self.struct:
            raise ValueError('need structure,')
        magmom = self.incar.get('magmom', [])
        lsorb = self.incar.get('lsorbit', False)
        lnonc = self.incar.get('lnoncollinear', False)
        ok = True
        nmag = len(magmom)
        nsit = len(self.struct.sites)
        if lsorb:
            if lnonc:
                if not nmag == 3*nsit:
                    ok = False
            else:
                if not nmag == nsit:
                    ok = False
        else:
            if not nmag == nsit:
                ok = False
        return ok

