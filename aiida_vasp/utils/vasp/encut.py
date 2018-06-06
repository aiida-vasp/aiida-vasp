"""Helper for creating useful ENCUT INCAR tags."""
import enum
from aiida.common.utils import classproperty

from aiida_vasp.io.potcar import MultiPotcarIo
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.utils.vasp.incar_param import AbstractIncarParam


# pylint: disable=too-few-public-methods
class EncutFlags(enum.Enum):
    MIN_ENMAX = 'min(Enmax)'
    MAX_ENMAX = 'max(Enmax)'


class Encut(AbstractIncarParam):
    """
    Helper for creating useful ENCUT INCAR tags.

    In most cases when the ENCUT tag should be given it should depend on the set of POTCARs that will be used.
    This class provides an interface to grab "factor * min(ENMAX)" or "factor * max(ENMAX)" using the inputs
    one needs for a calculation anyway.

    The set of potcars is determined by calling ``PotcarData.get_potcars_from_structure()``. from those,
    the maximum or minimum ENMAX is grabbed.

    :param structure: StructureData, passed to ``PotcarData.get_potcars_from_structure()``
    :param potcar_family: String, passed to ``PotcarData.get_potcars_from_structure()``
    :param potcar_mapping: Dictionary, passed to ``PotcarData.get_potcars_from_structure()``
    :param min_or_max: EncutFlags.MIN_ENMAX or EncutFlags.MAX_ENMAX
    :param factor: a factor with which the min or max ENMAX should be multiplied
    """

    def __init__(self, structure, potcar_family, potcar_mapping, min_or_max=EncutFlags.MAX_ENMAX, factor=1):
        if not isinstance(min_or_max, EncutFlags):
            raise ValueError('min_or_max parameter must be of type {}'.format(EncutFlags))
        self._structure = structure
        self._potcar_family = potcar_family
        self._multipotcar = MultiPotcarIo.from_structure(structure,
                                                         get_data_class('potcar').get_potcars_from_structure(
                                                             structure=structure, family_name=potcar_family, mapping=potcar_mapping))
        self._min_or_max = min_or_max
        if min_or_max == EncutFlags.MIN_ENMAX:
            self._base_encut = min([potcar.pymatgen.enmax for potcar in self._multipotcar.potcars])
        else:
            self._base_encut = max([potcar.pymatgen.enmax for potcar in self._multipotcar.potcars])
        self._factor = factor

    @classproperty
    def name(self):
        return 'ENCUT'

    @property
    def value(self):
        return self._base_encut * self._factor

    def clean(self, incar_dict):
        return [], []

    def validate(self):
        pass

    @property
    def info(self):
        msg = '{factor}x{minmax}, unit: eV'
        factor = '' if self._factor == 1 else '{}x'.format(self._factor)

        return msg.format(factor, self._min_or_max.value)
