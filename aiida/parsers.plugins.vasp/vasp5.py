from aiida.parsers.plugins.vasp.base import BaseParser
from aiida.tools.codespecific.vasp.io.eigenval import EigParser

class Vasp5Parser(BaseParser):
    '''
    Parses all Vasp 5 calculations.
    '''
    def parser_with_retrieved(self, retrieved):
        self.check_state()
        self.out_folder = self.get_folder(retrieved)
        if not self.out_folder:
            return self.result(success=False)
        outcar = self.get_file(OUTCAR)
        if not outcar:
            self.logger.error('OUTCAR not found, ' +
                    'look at the scheduler output for troubleshooting')
            return self.result(success=False)

        new_nodes = []
        bands, kpout = self.read_eigenval()
        new_nodes.append(('bands', bands))
        if kpout:
            new_nodes.append(('kpoints', kpout))

    def read_eigenval(self):
        incar = self._calc.incar.get_dict()
        nsw = incar.get('nsw', 0)
        ibrion = incar.get('ibrion', (nsw in [0, 1]) and -1 or 0)
        kpoints = self._calc.get_inputdata_dict().get('kpoints')
        eig = self.get_file('EIGENVALUE')
        header, kp, bs = EigParser.parse_eigenval(eig)
        bsnode = DataFactory('array.bands')()
        bsnode.set_bands(bs, occupations=None) # TODO: set occupations
        kpout = None
        structure = self._calc.structure
        if kpoints:
            bsnode.set_kpointsdata(kpoints)
        else:
            kpout = DataFactory('array.kpoints')()
            kpout.set_cell_from_structure(structure)
            kpout.set_kpoints(kp[:,:3], weights=kp[:,3])
            bsnode.set_kpointsdata(kpout)
        if ibrion in [-1, 1, 2]:
            bsnode.set_cell_from_structure(structure)
        else: # TODO: read contcar and associate with bsnode
            self.logger.info('dynamic vasp simulation detected ' +
                    '(IBRION = %s) -> EIGENVAL refers to CONTCAR, ' +
                    'not POSCAR, the bands node will not have the ' +
                    'cell attribute set.')
        return bsnode, kpout
