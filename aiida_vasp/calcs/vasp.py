from base import VaspCalcBase
from nscf import NscfCalculation
from aiida.orm import DataFactory
from aiida.common.utils import classproperty


class VaspCalculation(NscfCalculation):
    '''
    General-purpose VASP calculation, retrieves everything,
    so if storage space is a concern, consider deriving from it
    and overriding the retrieve_list in the subclass so only
    the necessary files are retrieved from the server.
    '''
    default_parser = 'vasp.vasp'

    def _prepare_for_submission(self, tempfolder, inputdict):
        '''retrieve all output files potentially created by VASP'''
        calcinfo = super(VaspCalculation, self)._prepare_for_submission(
            tempfolder, inputdict)
        calcinfo.retrieve_list = VaspCalcBase.max_retrieve_list()
        return calcinfo

    def verify_inputs(self, inputdict, *args, **kwargs):
        # ~ notset_msg = 'input not set: %s'
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

    @classmethod
    def _get_paw_linkname(cls, kind):
        return 'paw_%s' % kind

    @property
    def _parameters(self):
        return {k.lower(): v for
                k, v in self.inp.parameters.get_dict().iteritems()}

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
            (py:method::VaspCalculation.use_charge_densities),
            False otherwise
        needs 'parameters' input to be set
        (py:method::VaspCalculation.use_parameters)
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
            True if a wavecar file must be used
            (py:method::VaspCalculation.use_wavefunctions),
            False otherwise
        needs 'parameters' input to be set
        (py:method::VaspCalculation.use_parameters)
        '''
        istrt_d = self.get_inputs_dict().get('wavefunctions') and 1 or 0
        istart = self._parameters.get('istart', istrt_d)
        if istart in [1, 2, 3]:
            return True
        else:
            return False

    @classmethod
    def new_parameters(self, **kwargs):
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
        py:class:`~aiida_vasp.calcs.base.CalcMeta` ref CalcMeta
        pick up internal parameters from the class body
        and insert them
        '''
        super(VaspCalculation, self)._init_internal_params()
        self._update_internal_params()
