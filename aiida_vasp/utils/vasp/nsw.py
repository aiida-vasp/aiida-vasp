"""Utilities to help choose useful NSW values."""
from aiida.common.utils import classproperty

from aiida_vasp.utils.vasp.incar_param import AbstractIncarParam
from aiida_vasp.utils.vasp.ibrion import IbrionFlags, Ibrion


class Nsw(AbstractIncarParam):
    """
    Helper for using the NSW INCAR tag to get a value that makes sense considering related tags.

    :param value: the desired number of ionic steps
    :param ibrion: the IBRION tag value
    :param incar_dict: a dict with normalized (lower case) incar key-value pairs
    """
    DEFAULT_STEP_NUM = 20

    def __init__(self, value=None, ibrion=None, incar_dict=None):
        if incar_dict:
            self._ibrion = Ibrion.from_value(incar_dict.get('ibrion'))
        if ibrion:
            if not isinstance(ibrion, Ibrion):
                ibrion = Ibrion.from_value(ibrion)
            self._ibrion = ibrion
        self._value = value
        if not self._value and self.should_be_zero():
            self._value = 0
        self.validate()

    def should_be_zero(self):
        if self._ibrion:
            return bool(self._ibrion.value == IbrionFlags.NO_UPDATE.value)
        return False

    @classproperty
    def name(self):
        return 'NSW'

    @property
    def value(self):
        return self._value

    @property
    def info(self):
        msg = 'Number of ionic steps that will be performed at maximum'
        if self.should_be_zero():
            msg = 'IBRION={}, no ionic updates will be performed, so ionic steps would be a waste of time'.format(self._ibrion.value)
        return msg

    def validate(self):
        if not self._value and not self._ibrion:
            raise ValueError('Please pass either `value` or one of `ibrion` or `incar`')
        if self._value > 0 and self.should_be_zero():
            raise ValueError('NSW should be zero when ionic updates are not performed')

    def clean(self, incar_dict):
        errors = []
        warnings = []
        ibrion = Ibrion.from_value(incar_dict.get('ibrion', None))
        if Nsw(value=self.value, ibrion=ibrion).should_be_zero():
            if self.value >= 1:
                warnings.append('ionic updates are switched off (IBRION) but ionic steps will be performed (NSW)')
        elif self.value == 0:
            warnings.append('ionic updates are on (IBRION) but no ionic steps will be performed (NSW=0)')
        return errors, warnings

    @classmethod
    def from_value(cls, value):
        return cls(value=value)
