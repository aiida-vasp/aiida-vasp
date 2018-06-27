"""Utilities to decide on appropriate value for IBRION tag and catch clashes with other tags."""
import enum

from aiida.common.utils import classproperty
from aiida_vasp.utils.vasp.incar_param import AbstractIncarParam


class IbrionFlags(enum.Enum):
    """Enum with more descriptive names for the various possible IBRION values."""
    NO_UPDATE = -1
    MOLECULAR_DYNAMICS = 0
    IONIC_RELAXATION_RMM_DIIS = 1
    IONIC_RELAXATION_CG = 2
    IONIC_RELAXATION_DAMPED_MD = 3
    HESSIAN_FD_FULL = 5
    HESSIAN_FD_SYMMETRIZED = 6
    HESSIAN_PT_FULL = 7
    HESSIAN_PT_SYMMETRIZED = 8
    IMPROVED_DIMER_METHOD = 44

    @classproperty
    def ionic_relaxation_options(self):
        return [self.IONIC_RELAXATION_RMM_DIIS, self.IONIC_RELAXATION_CG, self.IONIC_RELAXATION_DAMPED_MD]


class Ibrion(AbstractIncarParam):
    """
    Set and validate IBRION INCAR tag.

    :param:ion_updates: one of the options in aiida_vasp.utils.vasp.ibrion.IbrionFlags, default NO_UPDATE, determines the actual value.

    Use an enum to set the value, to avoid typos and increase code clarity. Validates upon construction and provides a ``.clean()`` method,
    which allows scanning a dictionary of incar files for problems related to the current IBRION value.

    Usage::

        from aiida_vasp.utils.vasp.ibrion import IbrionFlags, Ibrion

        incar = {...}
        ibrion = Ibrion(ion_updates=IbrionFlags.MOLECULAR_DYNAMICS)
        incar.update(ibrion.param)
        errors, warnings = ibrion.clean(incar)
        for message in errors + warnings:
            print message  # or log, etc
        if errors:
            raise ValueError(errors[0])
    """

    IBRION_EXPLANATIONS = {
        IbrionFlags.NO_UPDATE: 'Ions are not moved',
        IbrionFlags.MOLECULAR_DYNAMICS: 'Standard ab-initio molecular dynamics',
        IbrionFlags.IONIC_RELAXATION_RMM_DIIS: 'Ionic relaxation using RMM-DIIS (quasi-newton algorithm)',
        IbrionFlags.IONIC_RELAXATION_CG: 'Ionic relaxation using conjugate-gradient algorithm',
        IbrionFlags.IONIC_RELAXATION_DAMPED_MD: 'Ionic relaxation using damped second order equation of motion',
        IbrionFlags.HESSIAN_FD_FULL: 'Calculate Hessian using finite-differences, full trial displacements, selective dynamics supported',
        IbrionFlags.HESSIAN_FD_SYMMETRIZED: 'Calculate Hessian using finite-differences, symmetrized trial displacements',
        IbrionFlags.HESSIAN_PT_FULL: 'Calculate Hessian using perturbation theory, full trial displacements',
        IbrionFlags.HESSIAN_PT_SYMMETRIZED: 'Calculate Hessian using perturbation theory, symmetrized trial displacements',
        IbrionFlags.IMPROVED_DIMER_METHOD: 'Move ions using the Improved Dimer Method'
    }

    def __init__(self, ion_updates=IbrionFlags.NO_UPDATE):
        if not isinstance(ion_updates, IbrionFlags):
            raise ValueError('ion_updates parameter must be of type {}'.format(IbrionFlags))
        self._ion_updates = ion_updates
        self.validate()

    @classproperty
    def name(self):
        return 'IBRION'

    @staticmethod
    def _validate_value(value):
        allowed_values = [flag.value for flag in IbrionFlags]
        if value not in allowed_values:
            raise ValueError('Value {} not supported for IBRION'.format(value))

    def validate(self):
        self._validate_value(self.value)

    @property
    def value(self):
        return self._ion_updates.value

    def clean(self, incar_dict):
        errors = []
        warnings = []
        if self._ion_updates != IbrionFlags.NO_UPDATE:
            if incar_dict.get('nsw', 0) == 0:
                warnings.append('NSW (# ionic steps) was not set or set to 0 when updating ions ({})'.format(self))
        if self._ion_updates == IbrionFlags.MOLECULAR_DYNAMICS:
            if 'potim' not in incar_dict:
                errors.append('POTIM must be given when running MD ({})'.format(self))
        if self._ion_updates in IbrionFlags.ionic_relaxation_options:  # pylint: disable=unsupported-membership-test
            if 'potim' not in incar_dict:
                warnings.append('POTIM was not given when running ionic relaxation ({})'.format(self))
        return errors, warnings

    @property
    def info(self):
        return self.IBRION_EXPLANATIONS[self._ion_updates]

    @classmethod
    def from_value(cls, value):
        cls._validate_value(value)
        ibrion_flag = [flag for flag in IbrionFlags if flag.value == value][0]
        return cls(ion_updates=ibrion_flag)
