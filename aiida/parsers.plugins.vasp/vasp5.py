from aiida.parsers.plugins.vasp.base import BaseParser
from aiida.tools.codespecific.vasp.io.eigenval import EigParser
from aiida.tools.codespecific.vasp.io.vasprun import VasprunParser
from aiida.tools.codespecific.vasp.io.doscar import DosParser
from aiida.orm import DataFactory


class Vasp5Parser(BaseParser):
    '''
    Parses all Vasp 5 calculations.
    '''
    def parser_with_retrieved(self, retrieved):
        self.check_state()
        self.out_folder = self.get_folder(retrieved)
        if not self.out_folder:
            return self.result(success=False)
        outcar = self.get_file('OUTCAR')
        if not outcar:
            self.logger.error(
                'OUTCAR not found, ' +
                'look at the scheduler output for troubleshooting')
            return self.result(success=False)

        self.vrp = self.read_run()

        bands, kpout = self.read_eigenval()

        structure = None  # get output structure if not static
        if self.vrp.is_md or self.vrp.is_relaxation:
            structure = self.read_cont()

        if self.vrp.is_md:  # set cell from input or output structure
            cellst = structure
        else:
            cellst = self._calc.inp.structure
        bands.set_cell_from_structure(cellst)
        kpout.set_cell_from_structure(cellst)

        self.dcp = self.read_dos()
        dosnode = self.get_dos_node(self.vrp, self.dcp)

        self.set_bands(bands)  # append output nodes
        if not self._calc.inp.kpoints:
            self.set_kpoints(kpout)
        if structure:
            self.set_structure(structure)

        if self.vrp.is_sc:  # add chgcar ouput node if selfconsistent run
            chgnode = self.get_chgcar()
            self.set_chgcar(chgnode)

        self.set_dos(dosnode)

        return self.new_nodes

    def read_run(self):
        '''Read vasprun.xml'''
        vasprun = self.get_file('vasprun.xml')
        if not vasprun:
            self.logger.warning('no vasprun.xml found')
            return None
        return VasprunParser(vasprun)

    def read_dos(self):
        '''read DOSCAR for more accurate tdos and pdos'''
        doscar = self.get_file('DOSCAR')
        if not doscar:
            self.logger.warning('no DOSCAR found')
            return None
        return DosParser(doscar)

    def get_dos_node(self, vrp, dcp):
        '''
        takes VasprunParser and DosParser objects
        and returns a doscar array node
        '''
        dosnode = DataFactory('array')()
        pdos = vrp.pdos.copy()
        for i, name in enumerate(vrp.pdos.dtype.names[1:]):
            ns = vrp.pdos.shape[1]
            # ~ pdos[name] = dcp[:, :, i+1:i+1+ns].transpose(0,2,1)
            cur = dcp.pdos[:, :, i+1:i+1+ns].transpose(0,2,1)
            cond = vrp.pdos[name] < 0.1
            pdos[name] = np.where(cond, cur, vrp.pdos[name])
        ns = 1
        if dcp.tdos.shape == 5:
            ns = 2
        tdos = vrp.tdos[:ns, :].copy()
        for i, name in enumerate(vrp.tdos.dtype.names[1:]):
            cur = dcp.tdos[:,i+1:i+1+ns].transpose()
            cond = vrp.tdos[name] < 0.1
            tdos[name] = np.where(cond, cur, vrp.tdos[name])
        dosnode.set_array('pdos', pdos)
        dosnode.set_array('tdos', tdos)
        return dosnode

    def read_cont(self):
        '''read CONTCAR for output structure'''
        from ase.io.vasp import read_vasp
        structure = DataFactory('structure')()
        cont = self.get_file('CONTCAR')
        if not cont:
            self.logger.info('CONTCAR not found!')
            return None
        structure.set_ase(read_vasp(cont))
        return structure

    def read_eigenval(self):
        '''
        Create a bands and a kpoints node from values in eigenvalue.

        returns: bsnode, kpout
        - bsnode: BandsData containing eigenvalues from EIGENVAL
                and occupations from vasprun.xml
        - kpout: KpointsData containing kpoints from EIGENVAL,

        both bsnode as well as kpnode come with cell unset
        '''
        eig = self.get_file('EIGENVAL')
        if not eig:
            self.logger.warning('EIGENVAL not found')
            return None, None
        header, kp, bs = EigParser.parse_eigenval(eig)
        bsnode = DataFactory('array.bands')()
        bsnode.set_bands(bs, occupations=self.vrp.occupations)
        kpout = DataFactory('array.kpoints')()
        kpout.set_kpoints(kp[:, :3], weights=kp[:, 3],
                          cartesian=header['cartesian'])
        return bsnode, kpout

    def get_chgcar(self):
        chgc = self.get_file('CHGCAR')
        chgnode = DataFactory('singlefile')
        chgnode.set_file(chgc)
        return chgnode

    def get_ouptput(self):
        output = DataFactory('parameter')()
        return output

    def set_bands(self, node):
        self.add_node('bands', node)

    def set_kpoints(self, node):
        self.add_node('kpoints', node)

    def set_chgcar(self, node):
        self.add_node('chgcar', node)

    def set_structure(self, node):
        self.add_node('structure', node)

    def set_dos(self, node):
        self.add_node('dos', node)
