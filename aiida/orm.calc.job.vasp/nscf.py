from base import VaspCalcBase, Input
from aiida.orm import DataFactory
from aiida.common.utils import classproperty


class NscfCalculation(VaspCalcBase):
    '''
    Calculation written and tested for vasp 5.3.5
    '''
    settings = Input(types='parameter',
                     doc='parameter node: parameters to be written to ' +
                     'the INCAR file')
    structure = Input(types=['structure', 'cif'],
                      doc='aiida structure node: will be converted to POSCAR')
    paw = Input(types='vasp.paw',
                doc='PAW nodes for each kind of element in the material\n' +
                'will be concatenated into POTCAR',
                param='kind')
    kpoints = Input(types='array.kpoints', doc='aiida kpoints node: ' +
                    'will be written to KPOINTS file')
    charge_density = Input(types='vasp.chargedensity',
                           doc='chargedensity node: should be obtained \n' +
                           'from the output of a selfconsistent ' +
                           'NscfCalculation (written to CHGCAR)')
    wavefunctions = Input(types='vasp.wavefun',
                          doc='wavefunction node: to speed up convergence ' +
                          'for continuation jobs')
    default_parser = 'vasp.nscf'

    def _prepare_for_submission(self, tempfolder, inputdict):
        '''changes the retrieve_list to retrieve only files needed
        to continue with nonscf runs'''
        calcinfo = super(NscfCalculation, self)._prepare_for_submission(
            tempfolder, inputdict)
        calcinfo.retrieve_list = [
            'EIGENVAL', 'DOSCAR', 'OUTCAR', 'vasprun.xml', 'wannier90*']
        return calcinfo

    def write_incar(self, inputdict, dst):
        '''
        converts from settings node (ParameterData) to INCAR format
        and writes to dst
        :param
            inputdict: required by baseclass
            dst: absolute path of the file to write to
        '''
        from incar import dict_to_incar
        with open(dst, 'w') as incar:
            incar.write(dict_to_incar(self.inp.settings.get_dict()))

    def write_poscar(self, inputdict, dst):
        '''
        converts from structures node (StructureData) to POSCAR format
        and writes to dst
        :param
            inputdict: required by baseclass
            dst: absolute path of the file to write to
        '''
        from ase.io.vasp import write_vasp
        with open(dst, 'w') as poscar:
            write_vasp(poscar, self.inp.structure.get_ase(), vasp5=True)

    def write_potcar(self, inputdict, dst):
        '''
        concatenatest multiple paw files into a POTCAR
        :param
            inputdict: required by baseclass
            dst: absolute path of the file to write to
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
        :param
            inputdict: required by baseclass
            dst: absolute path of the file to write to
        '''
        from base import kplitemp, kpltemp
        kp = self.inp.kpoints
        kpl, weights = kp.get_kpoints(also_weights=True)
        kw = zip(kpl, weights)
        with open(dst, 'w') as kpoints:
            kpls = '\n'.join(
                [kplitemp.format(k=k[0], w=k[1]) for k in kw])
            kps = kpltemp.format(N=len(kw), klist=kpls)
            kpoints.write(kps)

    def write_additional(self, tempfolder, inputdict):
        super(NscfCalculation, self).write_additional(
            tempfolder, inputdict)
        if self._need_chgd():
            chgcar = tempfolder.get_abs_path('CHGCAR')
            self.write_chgcar(inputdict, chgcar)
        if self._need_wfn():
            wavecar = tempfolder.get_abs_path('WAVECAR')
            self.write_wavecar(inputdict, wavecar)

    def write_chgcar(self, inputdict, dst):
        import shutil
        shutil.copyfile(self.inp.charge_density.get_file_abs_path(), dst)

    def write_wavecar(self, inputdict, dst):
        import shutil
        shutil.copyfile(self.inp.wavefunctions.get_file_abs_path(), dst)

    def verify_inputs(self, inputdict, *args, **kwargs):
        # ~ notset_msg = 'input not set: %s'
        super(NscfCalculation, self).verify_inputs(inputdict, *args, **kwargs)
        self.check_input(inputdict, 'settings')
        self.check_input(inputdict, 'structure')
        if 'elements' not in self.attrs():
            self._prestore()
        for kind in self.elements:
            self.check_input(inputdict, self._get_paw_linkname(kind))
        self.check_input(inputdict, 'kpoints', self._need_kp)
        self.check_input(inputdict, 'charge_density', self._need_chgd)
        self.check_input(inputdict, 'wavefunctions', self._need_wfn)
        kpoints = inputdict['kpoints']
        try:
            kp, w = kpoints.get_kpoints(also_weights=True)
        except AttributeError as e:
            print 'NscfCalculation needs a list of kpoints (not a mesh)!'
            raise e

    @classmethod
    def _get_paw_linkname(cls, kind):
        return 'paw_%s' % kind

    @property
    def _settings(self):
        return {k.lower(): v for
                k, v in self.inp.settings.get_dict().iteritems()}

    def _prestore(self):
        '''
        set attributes prior to storing
        '''
        super(NscfCalculation, self)._prestore()
        self._set_attr('input_kp_used', self._need_kp())
        self._set_attr('input_chgd_used', self._need_chgd())
        self._set_attr('input_wfn_used', self._need_wfn())
        self._set_attr('elements', list(set(
            self.inp.structure.get_ase().get_chemical_symbols())))

    def _need_kp(self):
        return True

    def _need_chgd(self):
        '''
        Test wether an charge_densities input is needed or not.
        :return output:
            True if a chgcar file must be used (py:method::NscfCalculation.use_charge_densities),
            False otherwise
        needs 'settings' input to be set (py:method::NscfCalculation.use_settings)
        '''
        ichrg_d = self._need_wfn() and 0 or 2
        icharg = self._settings.get('icharg', ichrg_d)
        if icharg in [1, 11]:
            return True
        else:
            return False

    def _need_wfn(self):
        '''
        Test wether a wavefunctions input is needed or not.
        :return output:
            True if a wavecar file must be used(py:method::NscfCalculation.use_wavefunctions),
            False otherwise
        needs 'settings' input to be set (py:method::NscfCalculation.use_settings)
        '''
        # ~ nsw = self._settings.get('nsw', 0)
        # ~ ibrion_d = nsw in [0, 1] and -1 or 0
        # ~ ibrion = self._settings.get('ibrion', ibrion_d)
        istrt_d = self.get_inputs_dict().get('wavefunctions') and 1 or 0
        istart = self._settings.get('istart', istrt_d)
        if istart in [1, 2, 3]:
            return True
        else:
            return False

    @classmethod
    def new_settings(self, **kwargs):
        return DataFactory('parameter')(**kwargs)

    @classmethod
    def new_structure(self, **kwargs):
        return DataFactory('structure')(**kwargs)

    @classmethod
    def new_kpoints(self, **kwargs):
        return DataFactory('array.kpoints')(**kwargs)

    @classmethod
    def new_charge_density(self, **kwargs):
        return DataFactory('vasp.chargedensity')(**kwargs)

    @classmethod
    def new_wavefunctions(self, **kwargs):
        return DataFactory('vasp.wavefun')(**kwargs)

    def new_wannier_settings(self, **kwargs):
        return DataFactory('parameter')(**kwargs)

    def new_wannier_data(self, **kwargs):
        return DataFactory('vasp.archive')(**kwargs)

    @classmethod
    def load_paw(self, *args, **kwargs):
        return self.Paw.load_paw(*args, **kwargs)[0]

    @classproperty
    def Paw(self):
        return DataFactory('vasp.paw')

    @property
    def input_kp_used(self):
        return self.get_attr('input_kp_used')

    @property
    def input_charge_density_used(self):
        return self.get_attr('input_chgd_used')

    @property
    def input_wavefunctions_used(self):
        return self.get_attr('input_wfn_used')

    @property
    def elements(self):
        return self.get_attr('elements')

    def _init_internal_params(self):
        '''
        let the metaclass py:class:`~aiida.orm.calculation.job.vasp.base.CalcMeta` ref CalcMeta
        pick up internal parameters from the class body and insert them
        '''
        super(NscfCalculation, self)._init_internal_params()
        self._update_internal_params()
