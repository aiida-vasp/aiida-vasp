from base import BasicCalculation, Input
from aiida.orm import DataFactory
from aiida.common.utils import classproperty


class NscfCalculation(BasicCalculation):
    '''
    Calculation written and tested for vasp 5.3.5
    '''
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
        calcinfo.retrieve_list.extend([ 'EIGENVAL', 'DOSCAR'])
        calcinfo.retrieve_list.extend(['wannier90.win',
                                       'wannier90.mmn',
                                       'wannier90.amn',
                                       'wannier90.eig'])
        return calcinfo

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
        let the metaclass py:class:`~aiida.orm.calculation.job.vasp.base.CalcMeta` ref CalcMeta
        pick up internal parameters from the class body and insert them
        '''
        super(NscfCalculation, self)._init_internal_params()
        self._update_internal_params()
