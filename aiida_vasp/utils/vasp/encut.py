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
    VALUE_ENTRY = '<given value>'


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
    :param strategy: EncutFlags.VALUE_ENTRY for manual entry of a value or
        EncutFlags.MIN_ENMAX, EncutFlags.MAX_ENMAX to automatically retrieve the value from a set of POTCARs
    :param factor: a factor with which the min or max ENMAX should be multiplied
    :param value: only in combination with strategy VALUE_ENTRY: manually enter a value
    """

    def __init__(self, strategy, **kwargs):
        if not isinstance(strategy, EncutFlags):
            raise ValueError('strategy parameter must be of type {}'.format(EncutFlags))
        self._strategy = strategy

        self._base_encut = kwargs.pop('value', None)
        self._factor = kwargs.pop('factor', 1)
        self._structure = kwargs.pop('structure', None)
        self._potcar_family = kwargs.pop('potcar_family', None)
        self._potcar_mapping = kwargs.pop('potcar_mapping', None)
        self._multipotcar = None
        self.validate()
        if self._strategy != EncutFlags.VALUE_ENTRY:
            self._init_from_structure(**kwargs)

    def _init_from_structure(self):
        """Get the value from the POTCARs depending on the `strategy` keyword passed in the ctor."""
        self._multipotcar = MultiPotcarIo.from_structure(
            self._structure,
            get_data_class('vasp.potcar').get_potcars_dict(
                elements=self._structure.get_kind_names(), family_name=self._potcar_family, mapping=self._potcar_mapping))
        if self._strategy == EncutFlags.MIN_ENMAX:
            self._base_encut = min([potcar.pymatgen.enmax for potcar in self._multipotcar.potcars])
        else:
            self._base_encut = max([potcar.pymatgen.enmax for potcar in self._multipotcar.potcars])

        if not self._base_encut:
            raise self.EncutInferringError('An unexpected error occurred while inferring ENCUT with strategy {}'.format(self._strategy))

    @classproperty
    def name(self):
        return 'ENCUT'

    @property
    def value(self):
        return self._base_encut * self._factor

    def clean(self, incar_dict):
        return [], []

    def validate(self):
        if self._strategy == EncutFlags.VALUE_ENTRY:
            if not self._base_encut:
                raise ValueError('missing kwarg "value" for strategy {}'.format(self._strategy))
        else:
            if not self._structure:
                raise ValueError('missing kwarg "structure" for strategy {}'.format(self._strategy))
            if not self._potcar_family:
                raise ValueError('missing kwarg "potcar_family" for strategy {}'.format(self._strategy))
            if not self._potcar_mapping:
                raise ValueError('missing kwarg "potcar_mapping" for strategy {}'.format(self._strategy))

    @property
    def info(self):
        msg = '{factor}x{strategy}, unit: eV'
        factor = '' if self._factor == 1 else '{}x'.format(self._factor)

        return msg.format(factor=factor, strategy=self._strategy.value)

    class EncutInferringError(Exception):
        """Raise when an error occurs while getting ENCUT from max or min(ENMAX) of the given potentials."""
        pass
