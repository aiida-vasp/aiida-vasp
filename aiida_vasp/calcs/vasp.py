try:
    from collections import ChainMap
except ImportError:
    from chainmap import ChainMap

from base import VaspCalcBase, Input
from aiida.orm import DataFactory


class VaspCalculation(VaspCalcBase):
    '''
    General-purpose VASP calculation. By default retrieves only the 'OUTCAR', 'vasprun.xml', 'EIGENVAL', 'DOSCAR' and Wannier90 input / output files, but additional retrieve files can be specified via the 'settings['ADDITIONAL_RETRIEVE_LIST']' input.
    '''
    default_parser = 'vasp.vasp'
    parameters = Input(types='parameter', doc='VASP INCAR parameters.')
    structure = Input(types=['structure', 'cif'])
    paw = Input(types='vasp.paw', param='kind')
    kpoints = Input(types='array.kpoints')
    settings = Input(types='parameter', doc='Additional settings for the calculation.')
    charge_density = Input(types='vasp.chargedensity', doc='chargedensity node: should be obtained from the output of a selfconsistent SCF calculation (written to CHGCAR)')
    wavefunctions = Input(types='vasp.wavefun', doc='wavefunction node: to speed up convergence for continuation jobs')

    _DEFAULT_PARAMETERS = {}
    _ALWAYS_RETRIEVE_LIST = ['OUTCAR', 'vasprun.xml', 'EIGENVAL', 'DOSCAR', ('wannier90*', '.', 0)]

    def _prepare_for_submission(self, tempfolder, inputdict):
        '''add EIGENVAL, DOSCAR, and all files starting with wannier90 to
        the list of files to be retrieved.'''
        calcinfo = super(VaspCalculation, self)._prepare_for_submission(
            tempfolder, inputdict)
        try:
            additional_retrieve_list = inputdict['settings'].get_attr('ADDITIONAL_RETRIEVE_LIST')
        except KeyError, AttributeError:
            additional_retrieve_list = []
        calcinfo.retrieve_list = list(set(
            self._ALWAYS_RETRIEVE_LIST + additional_retrieve_list
        ))
        return calcinfo

    def verify_inputs(self, inputdict, *args, **kwargs):
        super(VaspCalculation, self).verify_inputs(inputdict, *args, **kwargs)
        self.check_input(inputdict, 'parameters')
        self.check_input(inputdict, 'structure')
        if 'elements' not in self.attrs():
            self._prestore()
        for kind in self.elements:
            self.check_input(inputdict, self._get_paw_linkname(kind))
        self.check_input(inputdict, 'kpoints', self._need_kp)
        self.check_input(inputdict, 'charge_density', self._need_chgd)
        self.check_input(inputdict, 'wavefunctions', self._need_wfn)

    def _prestore(self):
        '''
        set attributes prior to storing
        '''
        super(VaspCalculation, self)._prestore()
        self._set_attr('elements', ordered_unique_list(
            self.inp.structure.get_ase().get_chemical_symbols()))

    @classmethod
    def _get_paw_linkname(cls, kind):
        '''required for storing multiple input paw nodes'''
        return 'paw_%s' % kind

    @property
    def _parameters(self):
        all_parameters = ChainMap(
            self.inp.parameters.get_dict(), self._DEFAULT_PARAMETERS
        )
        return {k.lower(): v for k, v in all_parameters.items()}

    @property
    def elements(self):
        return self.get_attr('elements')

    def _need_kp(self):
        '''
        return wether an input kpoints node is needed or not.
        :return output:
            True if input kpoints node is needed
            (py:method::VaspCalculation.use_kpoints),
            False otherwise
        needs 'parameters' input to be set
        (py:method::VaspCalculation.use_parameters)
        '''
        if 'kspacing' in self._parameters and 'kgamma' in self._parameters:
            return False
        else:
            return True

    def _need_chgd(self):
        '''
        Test wether an charge_densities input is needed or not.
        :return output:
            True if a chgcar file must be used
            (py:method::NscfCalculation.use_charge_densities),
            False otherwise
        needs 'parameters' input to be set
        (py:method::NscfCalculation.use_parameters)
        '''
        ichrg_d = self._need_wfn() and 0 or 2
        icharg = self._parameters.get('icharg', ichrg_d)
        if icharg in [1, 11]:
            return True
        else:
            return False

    def _need_wfn(self):
        '''
        Test wether a wavefunctions input is needed or not.
        :return output:
            True if a wavecar file must be
            used (py:method::NscfCalculation.use_wavefunctions),
            False otherwise
        needs 'parameters' input to be set
        (py:method::NscfCalculation.use_parameters)
        '''
        istrt_d = self.get_inputs_dict().get('wavefunctions') and 1 or 0
        istart = self._parameters.get('istart', istrt_d)
        if istart in [1, 2, 3]:
            return True
        else:
            return False

    def write_additional(self, tempfolder, inputdict):
        '''write CHGAR and WAVECAR files if needed'''
        super(VaspCalculation, self).write_additional(
            tempfolder, inputdict)
        if self._need_chgd():
            chgcar = tempfolder.get_abs_path('CHGCAR')
            self.write_chgcar(inputdict, chgcar)
        if self._need_wfn():
            wavecar = tempfolder.get_abs_path('WAVECAR')
            self.write_wavecar(inputdict, wavecar)

    def write_incar(self, inputdict, dst):
        '''
        Converts from parameters node (ParameterData) to INCAR format and writes to dst. Unless otherwise specified, the values specified in _DEFAULT_PARAMETERS are also written to the INCAR file.

        :param inputdict: required by baseclass
        :param dst: absolute path of the file to write to
        '''
        from ..utils.io.incar import dict_to_incar
        with open(dst, 'w') as incar:
            incar.write(dict_to_incar(
                ChainMap(self.inp.parameters.get_dict(), self._DEFAULT_PARAMETERS)
            ))

    def write_poscar(self, inputdict, dst):
        '''
        converts from structures node (StructureData) to POSCAR format
        and writes to dst

        :param inputdict: required by baseclass
        :param dst: absolute path of the file to write to
        '''
        from ase.io.vasp import write_vasp
        with open(dst, 'w') as poscar:
            write_vasp(poscar, self.inp.structure.get_ase(), vasp5=True)

    def write_potcar(self, inputdict, dst):
        '''
        Concatenates multiple paw files into a POTCAR

        :param inputdict: required by baseclass
        :param dst: absolute path of the file to write to
        '''
        import subprocess32 as sp
        catcom = ['cat']
        # ~ structure = inputdict['structure']
        # ~ structure = self.inp.structure
        # order the symbols according to order given in structure
        if 'elements' not in self.attrs():
            self._prestore()
        for kind in self.elements:
            paw = inputdict[self._get_paw_linkname(kind)]
            catcom.append(paw.get_abs_path('POTCAR'))
        # cat the pawdata nodes into the file
        with open(dst, 'w') as pc:
            sp.check_call(catcom, stdout=pc)

    def write_kpoints(self, inputdict, dst):
        '''
        converts from kpoints node (KpointsData) to KPOINTS format
        and writes to dst

        :param inputdict: required by baseclass
        :param dst: absolute path of the file to write to
        '''
        kp = self.inp.kpoints
        if kp.get_attrs().get('mesh'):
            self._write_kpoints_mesh(dst)
        elif kp.get_attrs().get('array|kpoints'):
            self._write_kpoints_list(dst)
        else:
            raise AttributeError('you supplied an empty kpoints node')

    def _write_kpoints_mesh(self, dst):
        kp = self.inp.kpoints
        mesh, offset = kp.get_kpoints_mesh()
        kpmtemp = (
            "Automatic mesh\n"
            "0\n"
            "Gamma\n"
            "{N[0]} {N[1]} {N[2]}\n"
            "{s[0]} {s[1]} {s[2]}\n"
        )
        with open(dst, 'w') as kpoints:
            kps = kpmtemp.format(N=mesh, s=offset)
            kpoints.write(kps)

    def _write_kpoints_list(self, dst):
        kp = self.inp.kpoints
        if 'array|weights' in kp.get_attrs():
            kpl, weights = kp.get_kpoints(also_weights=True)
        else:
            kpl = kp.get_kpoints()
            weights = [1.] * kpl.shape[0]
        kw = list(zip(kpl, weights))

        kpls = '\n'.join([
            '{k[0]} {k[1]} {k[2]} {w}'.format(k=k, w=w)
            for k, w in kw
        ])
        kps = (
            "Explicit list\n"
            "{N}\n"
            "Direct\n"
            "{klist}\n"
        ).format(N=len(kw), klist=kpls)
        with open(dst, 'w') as kpoints:
            kpoints.write(kps)

    def write_chgcar(self, inputdict, dst):
        import shutil
        shutil.copyfile(self.inp.charge_density.get_file_abs_path(), dst)

    def write_wavecar(self, inputdict, dst):
        import shutil
        shutil.copyfile(self.inp.wavefunctions.get_file_abs_path(), dst)

def ordered_unique_list(in_list):
    out_list = []
    for i in in_list:
        if i not in out_list:
            out_list.append(i)
    return out_list
