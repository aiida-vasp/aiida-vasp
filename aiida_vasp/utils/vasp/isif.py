"""Helper for using ISIF INCAR tag."""
import enum
from aiida.common.utils import classproperty

from aiida_vasp.utils.vasp.incar_param import AbstractIncarParam


# pylint: disable=too-few-public-methods
class IsifStressFlags(enum.Enum):
    NONE = 0
    FULL = 1
    TRACE_ONLY = 2


class Isif(AbstractIncarParam):
    """
    Helper for using ISIF INCAR tag, finds the correct value for a given set of degrees of freedom.

    :param calculate_stress: IsifStressFlags.Full | TRACE_ONLY | NONE
    :param vary_positions: bool, treat the ion positions as a degree of freedom?
    :param vary_cell_shape: bool, treat the cell shape as a degree of freedom?
    :param vary_cell_volume: bool, treat the cell volume as a degree of freedom?
    """

    ALLOWED_COMBINATIONS = {
        (IsifStressFlags.NONE, True, False, False): 0,
        (IsifStressFlags.TRACE_ONLY, True, False, False): 1,
        (IsifStressFlags.FULL, True, False, False): 2,
        (IsifStressFlags.FULL, True, True, True): 3,
        (IsifStressFlags.FULL, True, True, False): 4,
        (IsifStressFlags.FULL, False, True, False): 5,
        (IsifStressFlags.FULL, False, False, False): 7,
    }
    POS_KEY = 'Positions'
    SHAPE_KEY = 'Cell shape'
    VOL_KEY = 'Cell volume'

    def __init__(self, calculate_stress=IsifStressFlags.FULL, vary_positions=True, vary_cell_shape=False, vary_cell_volume=False):
        if not isinstance(calculate_stress, IsifStressFlags):
            raise ValueError('calculate_stress parameter must be of type {}'.format(IsifStressFlags))
        self._stress = calculate_stress
        self._dof = {self.POS_KEY: vary_positions, self.SHAPE_KEY: vary_cell_shape, self.VOL_KEY: vary_cell_volume}
        self.validate()

    @classproperty
    def name(self):
        return 'ISIF'

    @property
    def info(self):
        msg = 'Calculate stress tensor: {stress} | DoF: {dof} | fixed: {fixed}'
        stress = self._stress.name.capitalize()
        dof = ', '.join([key for key, value in self._dof.items() if value])
        fixed = ', '.join([key for key, value in self._dof.items() if not value])
        return msg.format(stress=stress, dof=dof, fixed=fixed)

    def _combination(self):
        return (self._stress, self._dof[self.POS_KEY], self._dof[self.SHAPE_KEY], self._dof[self.VOL_KEY])

    def validate(self):
        if self._combination() not in self.ALLOWED_COMBINATIONS:
            raise ValueError('Vasp does not allow the following combination of DoF: {}'.format(self.info))

    def clean(self, incar_dict):
        errors = []
        warnings = []
        if self._dof[self.VOL_KEY]:
            if 'encut' not in incar_dict:
                warnings.append(
                    'ENCUT was not given (should be set >= 1.3xmax(ENMAX)) when running cell volume relaxation ({})'.format(self))
        return errors, warnings

    @property
    def value(self):
        return self.ALLOWED_COMBINATIONS[self._combination()]

    @classmethod
    def from_value(cls, value):
        combination = [combo for combo, val in cls.ALLOWED_COMBINATIONS.items() if val == value][0]
        return cls(
            calculate_stress=combination[0], vary_positions=combination[1], vary_cell_shape=combination[2], vary_cell_volume=combination[2])
