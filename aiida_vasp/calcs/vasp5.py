# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""VASP - Calculation: Generic VASP5 calculation not relying on third party packages for preparing files"""
from aiida.orm import DataFactory
from aiida.common.utils import classproperty

from aiida_vasp.calcs.base import VaspCalcBase
from aiida_vasp.calcs.nscf import NscfCalculation
from aiida_vasp.calcs.wannier import WannierBase


class Vasp5Calculation(NscfCalculation, WannierBase):
    """
    General purpose VASP calculation, retrieves everything,
    so if storage space is a concern, consider deriving from it
    and overriding the retrieve_list in the subclass so only
    the necessary files are retrieved from the server.
    """
    default_parser = 'vasp.vasp5'

    def _prepare_for_submission(self, tempfolder, inputdict):
        """retrieve all output files potentially created by VASP"""
        calcinfo = super(Vasp5Calculation, self)._prepare_for_submission(
            tempfolder, inputdict)
        calcinfo.retrieve_list = VaspCalcBase.max_retrieve_list()
        return calcinfo

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
        return {
            k.lower(): v
            for k, v in self.inp.settings.get_dict().iteritems()
        }

    def _need_kp(self):
        """
        return wether an input kpoints node is needed or not.
        :return output:
            True if input kpoints node is needed
            (py:method::Vasp5Calculation.use_kpoints),
            False otherwise
        needs 'settings' input to be set
        (py:method::Vasp5Calculation.use_settings)
        """
        return bool('kspacing' in self._settings
                    and 'kgamma' in self._settings)

    def _need_chgd(self):
        """
        Test wether an charge_densities input is needed or not.
        :return output:
            True if a chgcar file must be used
            (py:method::Vasp5Calculation.use_charge_densities),
            False otherwise
        needs 'settings' input to be set
        (py:method::Vasp5Calculation.use_settings)
        """
        ichrg_d = 0 if self._need_wfn() else 2
        icharg = self._settings.get('icharg', ichrg_d)
        return bool(icharg in [1, 11])

    def _need_wfn(self):
        """
        Test wether a wavefunctions input is needed or not.
        :return output:
            True if a wavecar file must be used
            (py:method::Vasp5Calculation.use_wavefunctions),
            False otherwise
        needs 'settings' input to be set
        (py:method::Vasp5Calculation.use_settings)
        """
        istrt_d = 1 if self.get_inputs_dict().get('wavefunctions') else 0
        istart = self._settings.get('istart', istrt_d)
        return bool(istart in [1, 2, 3])

    @classmethod
    def new_settings(cls, **kwargs):
        return DataFactory('parameter')(**kwargs)

    @classmethod
    def new_structure(cls, **kwargs):
        return DataFactory('structure')(**kwargs)

    @classmethod
    def new_kpoints(cls, **kwargs):
        return DataFactory('array.kpoints')(**kwargs)

    @classmethod
    def new_charge_density(cls, **kwargs):
        return DataFactory('vasp.chargedensity')(**kwargs)

    @classmethod
    def new_wavefunctions(cls, **kwargs):
        return DataFactory('vasp.wavefun')(**kwargs)

    @classmethod
    def load_paw(cls, *args, **kwargs):
        return cls.PAW_CLS.load_paw(*args, **kwargs)[0]

    @classproperty
    def paw_cls(self):  # pylint: disable=no-self-use
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
        """
        let the metaclass
        py:class:`~aiida_vasp.calcs.base.CalcMeta` ref CalcMeta
        pick up internal parameters from the class body
        and insert them
        """
        super(Vasp5Calculation, self)._init_internal_params()
        self._update_internal_params()
