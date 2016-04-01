from nscf import NscfCalculation
from wannier import WannierBase
from aiida.orm import DataFactory
from aiida.common.utils import classproperty


class Vasp5Calculation(NscfCalculation, WannierBase):
    '''
    Calculation written and tested for vasp 5.3.5
    '''
    default_parser = 'vasp.vasp5'

    def write_additional(self, tempfolder, inputdict):
        super(Vasp5Calculation, self).write_additional(
            tempfolder, inputdict)
        if self._need_chgd():
            chgcar = tempfolder.get_abs_path('CHGCAR')
            self.write_chgcar(inputdict, chgcar)
        if self._need_wfn():
            wavecar = tempfolder.get_abs_path('WAVECAR')
            self.write_wavecar(inputdict, wavecar)

    def verify_inputs(self, inputdict, *args, **kwargs):
        # ~ notset_msg = 'input not set: %s'
        super(Vasp5Calculation, self).verify_inputs(inputdict, *args, **kwargs)
        self.check_input(inputdict, 'settings')
        self.check_input(inputdict, 'structure')
        if 'elements' not in self.attrs():
            self._prestore()
        for kind in self.elements:
            self.check_input(inputdict, self._get_paw_linkname(kind))
        self.check_input(inputdict, 'kpoints', self._need_kp)
        self.check_input(inputdict, 'charge_density', self._need_chgd)
        self.check_input(inputdict, 'wavefunctions', self._need_wfn)

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
        super(Vasp5Calculation, self)._prestore()

    def _need_kp(self):
        '''
        return wether an input kpoints node is needed or not.
        :return output:
            True if input kpoints node is needed
            (py:method::Vasp5Calculation.use_kpoints),
            False otherwise
        needs 'settings' input to be set
        (py:method::Vasp5Calculation.use_settings)
        '''
        if 'kspacing' in self._settings and 'kgamma' in self._settings:
            return False
        else:
            return True

    def _need_chgd(self):
        '''
        Test wether an charge_densities input is needed or not.
        :return output:
            True if a chgcar file must be used
            (py:method::Vasp5Calculation.use_charge_densities),
            False otherwise
        needs 'settings' input to be set
        (py:method::Vasp5Calculation.use_settings)
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
            True if a wavecar file must be used
            (py:method::Vasp5Calculation.use_wavefunctions),
            False otherwise
        needs 'settings' input to be set
        (py:method::Vasp5Calculation.use_settings)
        '''
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
        let the metaclass
        py:class:`~aiida.orm.calculation.job.vasp.base.CalcMeta` ref CalcMeta
        pick up internal parameters from the class body
        and insert them
        '''
        super(Vasp5Calculation, self)._init_internal_params()
        self._update_internal_params()
