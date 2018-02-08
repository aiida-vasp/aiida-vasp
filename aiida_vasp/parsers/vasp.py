#encoding: utf-8
"""AiiDA Parser for a aiida_vasp.VaspCalculation"""
import numpy as np
from aiida.orm import DataFactory

from aiida_vasp.io.doscar import DosParser
from aiida_vasp.io.eigenval import EigParser
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.vasprun import VasprunParser
from aiida_vasp.parsers.base import BaseParser


class VaspParser(BaseParser):
    """
    Parses all Vasp calculations.
    """

    def __init__(self, calc):
        super(VaspParser, self).__init__(calc)
        self.out_folder = None
        self.vrp = None
        self.dcp = None

    def parse_with_retrieved(self, retrieved):
        self.check_state()
        self.out_folder = self.get_folder(retrieved)
        if not self.out_folder:
            return self.result(success=False)
        outcar = self.get_file('OUTCAR')
        if not outcar:
            self.logger.error('OUTCAR not found, ' + 'look at the scheduler output for troubleshooting')
            return self.result(success=False)

        self.vrp = self.read_run()

        bands, kpout, structure = self.read_eigenval()

        self.dcp = self.read_dos()
        dosnode = self.get_dos_node(self.vrp, self.dcp)

        if bands:
            self.set_bands(bands)  # append output nodes

        if not kpout:
            kpout = self.read_ibzkpt()

        if kpout:
            self.set_kpoints(kpout)

        if structure:
            self.set_structure(structure)

        if self.vrp:
            if self.vrp.is_sc:  # add chgcar ouput node if selfconsistent run
                chgnode = self.get_chgcar()
                self.set_chgcar(chgnode)

            if self.vrp.is_sc:
                self.set_wavecar(self.get_wavecar())

        self.add_node('results', self.get_output())

        if dosnode:
            self.set_dos(dosnode)

        return self.result(success=True)

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

    @staticmethod
    def get_dos_node(vrp, dcp):
        """
        takes VasprunParser and DosParser objects
        and returns a doscar array node
        """
        if not vrp or not dcp:
            return None
        dosnode = DataFactory('array')()
        # vrp.pdos is a numpy array, and thus not directly bool-convertible
        if vrp.pdos.size > 0:
            pdos = vrp.pdos.copy()
            for i, name in enumerate(vrp.pdos.dtype.names[1:]):
                num_spins = vrp.pdos.shape[1]
                # ~ pdos[name] = dcp[:, :, i+1:i+1+ns].transpose(0,2,1)
                cur = dcp.pdos[:, :, i + 1:i + 1 + num_spins].transpose(0, 2, 1)
                cond = vrp.pdos[name] < 0.1
                pdos[name] = np.where(cond, cur, vrp.pdos[name])
            dosnode.set_array('pdos', pdos)
        num_spins = 1
        if dcp.tdos.shape[1] == 5:
            num_spins = 2
        tdos = vrp.tdos[:num_spins, :].copy()
        for i, name in enumerate(vrp.tdos.dtype.names[1:]):
            cur = dcp.tdos[:, i + 1:i + 1 + num_spins].transpose()
            cond = vrp.tdos[:num_spins, :][name] < 0.1
            tdos[name] = np.where(cond, cur, vrp.tdos[:num_spins, :][name])
        dosnode.set_array('tdos', tdos)
        return dosnode

    def read_cont(self):
        '''read CONTCAR for output structure'''
        from ase.io import read
        structure = DataFactory('structure')()
        cont = self.get_file('CONTCAR')
        if not cont:
            self.logger.info('CONTCAR not found!')
            return None
        structure.set_ase(read(cont, format='vasp'))
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
            return None, None, None
        _, kpoints, bands = EigParser.parse_eigenval(eig)
        bsnode = DataFactory('array.bands')()
        kpout = DataFactory('array.kpoints')()

        structure = None  # get output structure if not static
        if self.vrp.is_md or self.vrp.is_relaxation:
            structure = self.read_cont()

        if self.vrp.is_md:  # set cell from input or output structure
            cellst = structure
        else:
            cellst = self._calc.inp.structure
        bsnode.set_cell(cellst.get_ase().get_cell())
        kpout.set_cell(cellst.get_ase().get_cell())

        if self._calc.inp.kpoints.get_attrs().get('array|kpoints'):
            bsnode.set_kpointsdata(self._calc.inp.kpoints)
        if self._calc.inp.kpoints.labels:
            bsnode.labels = self._calc.inp.kpoints.labels
        else:
            bsnode.set_kpoints(kpoints[:, :3], weights=kpoints[:, 3], cartesian=False)
        bsnode.set_bands(bands, occupations=self.vrp.occupations)
        kpout.set_kpoints(kpoints[:, :3], weights=kpoints[:, 3], cartesian=False)
        return bsnode, kpout, structure

    def read_ibzkpt(self):
        """Create a DB Node for the IBZKPT file"""
        ibz = self.get_file('IBZKPT')
        if not ibz:
            self.logger.warning('IBZKPT not found')
            return None
        kpp = KpParser(ibz)
        kpout = DataFactory('array.kpoints')()
        kpout.set_kpoints(kpp.kpoints, weights=kpp.weights, cartesian=kpp.cartesian)
        return kpout

    def get_chgcar(self):
        """Create a DB Node for the CHGCAR file"""
        chgc = self.get_file('CHGCAR')
        if chgc is None:
            return None
        chgnode = DataFactory('vasp.chargedensity')()
        chgnode.set_file(chgc)
        return chgnode

    def get_wavecar(self):
        """Create a DB Node for the WAVECAR file"""
        wfn = self.get_file('WAVECAR')
        if wfn is None:
            return None
        wfnode = DataFactory('vasp.wavefun')()
        wfnode.set_file(wfn)
        return wfnode

    def get_output(self):
        output = DataFactory('parameter')()
        output.update_dict({'efermi': self.vrp.efermi})
        return output

    def set_bands(self, node):
        self.add_node('bands', node)

    def set_kpoints(self, node):
        self.add_node('kpoints', node)

    def set_chgcar(self, node):
        if node is not None:
            self.add_node('charge_density', node)

    def set_wavecar(self, node):
        if node is not None:
            self.add_node('wavefunctions', node)

    def set_structure(self, node):
        self.add_node('structure', node)

    def set_dos(self, node):
        self.add_node('dos', node)
