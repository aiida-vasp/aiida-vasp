"""
Parameter related utils

-----------------------
Contains utils and definitions that are used together with the parameters.
"""
import enum

from aiida.common.extendeddicts import AttributeDict
from aiida.plugins import DataFactory


def find_key_in_dicts(dictionary, supplied_key):
    """Find a key in a nested dictionary."""
    for key, value in dictionary.items():
        if key == supplied_key:
            yield value
        elif isinstance(value, dict):
            for result in find_key_in_dicts(value, supplied_key):
                yield result


class RelaxAlgoEnum(enum.IntEnum):
    """Encode values for algorithm descriptively in enum."""
    NO_UPDATE = -1
    IONIC_RELAXATION_RMM_DIIS = 1
    IONIC_RELAXATION_CG = 2


class RelaxModeEnum(enum.IntEnum):
    """
    Encode values for mode of relaxation descriptively in enum.

    Values can be found here: https://cms.mpi.univie.ac.at/wiki/index.php/ISIF
    """

    POS_ONLY = 2
    POS_SHAPE_VOL = 3
    POS_SHAPE = 4
    SHAPE_ONLY = 5
    SHAPE_VOL = 6
    VOL_ONLY = 7

    @classmethod
    def get_from_dof(cls, **kwargs):
        """Get the correct mode of relaxation for the given degrees of freedom."""
        RELAX_POSSIBILITIES = ('positions', 'shape', 'volume')  # pylint: disable=invalid-name
        dof = tuple(kwargs[i] for i in RELAX_POSSIBILITIES)
        value_from_dof = {
            (True, False, False): cls.POS_ONLY,
            (True, True, True): cls.POS_SHAPE_VOL,
            (True, True, False): cls.POS_SHAPE,
            (False, True, False): cls.SHAPE_ONLY,
            (False, True, True): cls.SHAPE_VOL,
            (False, False, True): cls.VOL_ONLY
        }
        try:
            return value_from_dof[dof]
        except KeyError:
            raise ValueError('Invalid combination for degrees of freedom: {}'.format(dict(zip(RELAX_POSSIBILITIES, dof))))


class ParametersMassage():
    """
    A class that contains all relevant massaging of the input parameters for VASP.

    The idea is that this class accepts the set input parameters from AiiDA (non code specifics), checks if any code specific
    parameters supplied are valid VASP input parameters (only rudimentary at this point, should also cross check and check types)
    and convert the AiiDA input parameters to VASP specific parameters. A set function needs to be developed for each parameter.
    This set function takes the AiiDA input and converts it. The parameter property should return ready to go parameters
    that can be dumped using the parsers in the respective CalcJob plugins.
    """

    def __init__(self, workchain, parameters):
        # First of all make sure parameters is a not a AiiDA Dict datatype
        self.exit_code = None
        self._workchain = workchain
        self._massage = AttributeDict()
        if isinstance(parameters, DataFactory('dict')):
            self._parameters = AttributeDict(parameters.get_dict())
        elif isinstance(parameters, AttributeDict):
            self._parameters = parameters
        else:
            raise TypeError('The supplied type: {} of parameters is not supported. '
                            'Supply either a Dict or an AttributeDict'.format(type(parameters)))
        self._load_valid_params()
        self._functions = ParameterSetFunctions(self._workchain, self._parameters, self._massage)
        self._prepare_parameters()
        self._check_parameters()

    def _load_valid_params(self):
        """Import a list of valid parameters for VASP. This is generated from the manual."""
        from os import path  # pylint: disable=import-outside-toplevel
        from yaml import safe_load  # pylint: disable=import-outside-toplevel
        with open(path.join(path.dirname(path.realpath(__file__)), 'tags.yml'), 'r') as file_handler:
            tags_data = safe_load(file_handler)
        self._valid_parameters = list(tags_data.keys())

    def _prepare_parameters(self):
        """Iterate over the valid parameters and call the set function associated with that parameter."""
        for key in self._valid_parameters:
            self._set(key)
            # We check after each parameter set if there is an exit code set on the WorkChain, if so, return
            if self.exit_code is not None:
                return

    def _check_parameters(self):
        """Make sure all the massaged values are to VASP spec."""
        if list(self._massage.keys()).sort() != self._valid_parameters.sort() and self.exit_code is None:
            self.exit_code = self._workchain.exit_codes.ERROR_INVALID_PARAMETER_DETECTED

    def _set(self, key):
        """Call the necessary function to set each parameter."""
        try:
            exit_code = getattr(self._functions, 'set_' + key)()
            if exit_code is not None:
                self.exit_code = exit_code
        except AttributeError:
            pass
        # If we find any raw code input key directly on parameter root, override whatever we have set until now
        # Also, make sure it is lowercase
        if self._parameters.get(key) or self._parameters.get(key.upper()):
            self._massage[key] = self._parameters[key]

    @property
    def parameters(self):
        """Return the massaged parameter set ready to go in VASP format."""
        return self._massage


class ParameterSetFunctions():
    """Container for the set functions that converts an AiiDA parameters to a code specific one."""

    def __init__(self, workchain, parameters, massage):
        self._parameters = parameters
        self._workchain = workchain
        self._massage = massage

    def set_ibrion(self):
        """Set which algorithm to use for ionic movements."""
        if self._relax():
            try:
                if self._parameters.relax.algo == 'cg':
                    self._massage.ibrion = RelaxAlgoEnum.IONIC_RELAXATION_CG.value
                elif self._parameters.relax.algo == 'rd':
                    self._massage.ibrion = RelaxAlgoEnum.IONIC_RELAXATION_RMM_DIIS.value
                else:
                    self._workchain.report('Invalid algo parameter: {}'.format(self._parameters.relax.algo))
                    return self._workchain.exit_codes.ERROR_INVALID_PARAMETER_DETECTED
            except AttributeError:
                self._workchain.report('Missing parameter: algo')
                return self._workchain.exit_codes.ERROR_MISSING_PARAMETER_DETECTED

        return None

    def set_ediffg(self):
        """Set the cutoff to use for relaxation."""
        if not self._relax():
            return
        energy_cutoff = False
        try:
            self._massage.ediffg = self._parameters.relax.energy_cutoff
            energy_cutoff = True
        except AttributeError:
            pass
        try:
            self._massage.ediffg = -abs(self._parameters.relax.force_cutoff)
            if energy_cutoff:
                self._workchain.report('User supplied both a force and an energy cutoff for the relaxation. Utilizing the force cutoff.')
        except AttributeError:
            pass

    def set_nsw(self):
        """Set the number of ionic steps to perform."""
        if self._relax():
            self._set_simple('nsw', self._parameters.relax.steps)

    def set_isif(self):
        """Set relaxation mode according to the chosen degrees of freedom."""
        positions = self._parameters.get('relax', {}).get('positions', False)
        shape = self._parameters.get('relax', {}).get('shape', False)
        volume = self._parameters.get('relax', {}).get('volume', False)
        if positions or shape or volume:
            self._massage.isif = RelaxModeEnum.get_from_dof(positions=positions, shape=shape, volume=volume).value

    def _relax(self):
        """Check if we have enabled relaxation."""
        return self._parameters.get('relax', {}).get('positions') or \
            self._parameters.get('relax', {}).get('shape') or \
            self._parameters.get('relax', {}).get('volume')

    def _set_simple(self, target, value):
        """Set basic parameter."""
        try:
            self._massage[target] = value
        except AttributeError:
            pass
