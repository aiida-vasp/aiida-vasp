from base import BasicCalculation, Input
from aiida.orm import DataFactory


class NscfCalculation(BasicCalculation):
    '''
    Runs VASP with precalculated (from scf run) wave functions and charge densities.
    Used to obtain bandstructures, DOS and wannier90 input files.
    '''
    charge_density = Input(types='vasp.chargedensity',
                           doc='chargedensity node: should be obtained \n' +
                           'from the output of a selfconsistent ' +
                           'NscfCalculation (written to CHGCAR)')
    wavefunctions = Input(types='vasp.wavefun',
                          doc='wavefunction node: to speed up convergence ' +
                          'for continuation jobs')
    default_parser = 'vasp.nscf'

    def verify_inputs(self, inputdict):
        super(NscfCalculation, self).verify_inputs(inputdict)
        self.check_input(inputdict, 'charge_density', self._need_chgd)
        self.check_input(inputdict, 'wavefunctions', self._need_wfn)
        if not self._need_chgd() and inputdict.get('charge_density'):
            msg = 'charge_density node given but '
            msg += '"icharg" key in settings not set '
            msg += 'to either 1 o 11. charge_density node not used --> .'
            msg += 'CHGCAR not written'
            self.logger.warning(msg)

    def _prepare_for_submission(self, tempfolder, inputdict):
        '''add EIGENVAL, DOSCAR, and all files starting with wannier90 to
        the list of files to be retrieved.'''
        calcinfo = super(NscfCalculation, self)._prepare_for_submission(
            tempfolder, inputdict)
        calcinfo.retrieve_list.extend(['EIGENVAL', 'DOSCAR'])
        calcinfo.retrieve_list.append(['wannier90*', '.', 0])
        return calcinfo

    def write_additional(self, tempfolder, inputdict):
        '''write CHGAR and WAVECAR files if needed'''
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

    def _prestore(self):
        '''
        set attributes prior to storing
        '''
        super(NscfCalculation, self)._prestore()
        self._set_attr('input_kp_used', self._need_kp())
        self._set_attr('input_chgd_used', self._need_chgd())
        self._set_attr('input_wfn_used', self._need_wfn())

    def _need_kp(self):
        return True

    def _need_chgd(self):
        '''
        Test wether an charge_densities input is needed or not.
        :return output:
            True if a chgcar file must be used
            (py:method::NscfCalculation.use_charge_densities),
            False otherwise
        needs 'settings' input to be set
        (py:method::NscfCalculation.use_settings)
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
            True if a wavecar file must be
            used (py:method::NscfCalculation.use_wavefunctions),
            False otherwise
        needs 'settings' input to be set
        (py:method::NscfCalculation.use_settings)
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
    def new_charge_density(self, **kwargs):
        return DataFactory('vasp.chargedensity')(**kwargs)

    @classmethod
    def new_wavefunctions(self, **kwargs):
        return DataFactory('vasp.wavefun')(**kwargs)

    def new_wannier_settings(self, **kwargs):
        return DataFactory('parameter')(**kwargs)

    def new_wannier_data(self, **kwargs):
        return DataFactory('vasp.archive')(**kwargs)

    @property
    def input_kp_used(self):
        return self.get_attr('input_kp_used')

    @property
    def input_charge_density_used(self):
        return self.get_attr('input_chgd_used')

    @property
    def input_wavefunctions_used(self):
        return self.get_attr('input_wfn_used')

    def _init_internal_params(self):
        '''
        let the metaclass
        py:class:`~aiida.orm.calculation.job.vasp.base.CalcMeta`
        ref CalcMeta
        pick up internal parameters from the class body and insert them
        '''
        super(NscfCalculation, self)._init_internal_params()
        self._update_internal_params()
