"""
Parameter related utils

-----------------------
Contains utils and definitions that are used together with the parameters.
"""

import enum

from aiida.common.extendeddicts import AttributeDict
from aiida.plugins import DataFactory
from aiida_vasp.utils.extended_dicts import update_nested_dict

FUNCTIONAL_PARAMETERS = {
    'charge': {
        'wave': True,
        'charge': True,
        'potential': True,
        'constant_charge': True,
        'constant_atomic': True
    },
    'smearing': {
        'mp': True,
        'gaussian': True,
        'fermi': True,
        'partial': True,
        'tetra': True
    }
}


class ChargeEnum(enum.IntEnum):
    """
    Encode values for the initial charge density.

    See: https://www.vasp.at/wiki/index.php/ICHARG
    """
    WAVE = 0
    CHARGE = 1
    ATOMIC = 2
    POTENTIAL = 4
    CONSTANT_CHARGE = 11
    CONSTANT_ATOMIC = 12


class IntSmearingEnum(enum.IntEnum):
    """
    Encode values for the smearing used during integration in reciprocal space.

    See: https://www.vasp.at/wiki/index.php/ISMEAR.
    """
    MP = 1  # pylint: disable=invalid-name
    GAUSSIAN = 0
    FERMI = -1
    PARTIAL = -2
    TETRA = -5


class OrbitEnum(enum.IntEnum):
    """
    Encode values for the projector information.

    See: https://www.vasp.at/wiki/index.php/LORBIT
    """
    ATOM = 0
    ATOM_LM = 1
    ATOM_LM_PHASE = 2
    NO_RWIGS_ATOM = 10
    NO_RWIGS_ATOM_LM = 11
    NO_RWIGS_ATOM_LM_PHASE = 12
    NO_RWIGS_ATOM_LM_PHASE_AUTO = 14
    ATOM_LM_WAVE = 5

    @classmethod
    def get_lorbit_from_combination(cls, **kwargs):
        """Get the correct mode of the projectors/decomposition."""
        combination = tuple(kwargs[i] for i in ['lm', 'phase', 'wigner_seitz_radius'])
        value_from_combinations = {
            (False, False, True): cls.ATOM,
            (True, False, True): cls.ATOM_LM,
            (True, True, True): cls.ATOM_LM_PHASE,
            (False, False, False): cls.NO_RWIGS_ATOM,
            (True, False, False): cls.NO_RWIGS_ATOM_LM,
            (True, True, False): cls.NO_RWIGS_ATOM_LM_PHASE,
            # Not supported, so also calculate lm decomposed
            (False, True, True): cls.ATOM_LM_PHASE,
            # Not supported, so also calculate lm decomposed
            (False, True, False): cls.NO_RWIGS_ATOM_LM_PHASE
        }
        return value_from_combinations[combination]


class RelaxAlgoEnum(enum.IntEnum):
    """
    Encode values for algorithm descriptively.

    See: https://www.vasp.at/wiki/index.php/ALGO
    """
    NO_UPDATE = -1
    IONIC_RELAXATION_RMM_DIIS = 1
    IONIC_RELAXATION_CG = 2


class RelaxModeEnum(enum.IntEnum):
    """
    Encode values for degrees of freedom mode of relaxation descriptively.

    See: https://cms.mpi.univie.ac.at/wiki/index.php/ISIF
    """

    POS_ONLY = 2
    POS_SHAPE_VOL = 3
    POS_SHAPE = 4
    SHAPE_ONLY = 5
    SHAPE_VOL = 6
    VOL_ONLY = 7

    @classmethod
    def get_isif_from_dof(cls, **kwargs):
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
        self._check_parameters()
        self._functions = ParameterSetFunctions(self._workchain, self._parameters, self._massage)
        self._set_parameters()
        self._set_vasp_parameters()
        self._validate_parameters()

    def _check_parameters(self):
        """Check that any directly supplied VASP parameters is valid."""

    def _load_valid_params(self):
        """Import a list of valid parameters for VASP. This is generated from the manual."""
        from os import path  # pylint: disable=import-outside-toplevel
        from yaml import safe_load  # pylint: disable=import-outside-toplevel
        with open(path.join(path.dirname(path.realpath(__file__)), 'tags.yml'), 'r') as file_handler:
            tags_data = safe_load(file_handler)
        self._valid_parameters = list(tags_data.keys())

    def _set_parameters(self):
        """Iterate over the valid parameters and call the set function associated with that parameter."""
        for key in self._valid_parameters:
            self._set(key)
            # We check after each parameter set if there is an exit code set on the WorkChain, if so, return
            if self.exit_code is not None:
                return

    def _set_vasp_parameters(self):
        """Set the any supplied override parameters."""
        try:
            for key, item in self._parameters.vasp.items():
                # Make sure key is lowercase on storage
                self._massage[key.lower()] = item
        except AttributeError:
            # Might not be supplied
            pass

    def _validate_parameters(self):
        """Make sure all the massaged values are to VASP spec."""
        if list(self._massage.keys()).sort() != self._valid_parameters.sort() and self.exit_code is None:
            self.exit_code = self._workchain.exit_codes.ERROR_INVALID_PARAMETER_DETECTED

    def _validate_vasp_parameters(self):
        """Validate the supplied VASP override parameters."""

    def _set(self, key):
        """Call the necessary function to set each parameter."""
        try:
            exit_code = getattr(self._functions, 'set_' + key)()
            if exit_code is not None:
                self.exit_code = exit_code
        except AttributeError:
            pass
        # If we find any raw code input key directly on parameter root, override whatever we have set until now
        # Note that the key may be in upper case, so we test both
        if key in self._parameters:
            self._massage[key] = self._parameters[key]
        elif key.upper() in self._parameters:
            self._massage[key] = self._parameters[key.upper()]

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
        """
        Set which algorithm to use for ionic movements.

        See: https://www.vasp.at/wiki/index.php/IBRION
        """
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
        """
        Set the cutoff to use for relaxation.

        See: https://www.vasp.at/wiki/index.php/EDIFFG
        """
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
        """
        Set the number of ionic steps to perform.

        See: https://www.vasp.at/wiki/index.php/NSW
        """
        if self._relax():
            self._set_simple('nsw', self._parameters.relax.steps)

    def set_isif(self):
        """
        Set relaxation mode according to the chosen degrees of freedom.

        See: https://www.vasp.at/wiki/index.php/ISIF
        """
        positions = self._parameters.get('relax', {}).get('positions', False)
        shape = self._parameters.get('relax', {}).get('shape', False)
        volume = self._parameters.get('relax', {}).get('volume', False)
        if positions or shape or volume:
            self._massage.isif = RelaxModeEnum.get_isif_from_dof(positions=positions, shape=shape, volume=volume).value

    def set_ismear(self):
        """
        Make sure we do not supply invalid integration methods when running explicit k-point grids.

        See: https://www.vasp.at/wiki/index.php/ISMEAR
        """

        try:
            if self._parameters.smearing.gaussian:
                self._set_simple('ismear', IntSmearingEnum.GAUSSIAN.value)
        except AttributeError:
            pass
        try:
            if self._parameters.smearing.fermi:
                self._set_simple('ismear', IntSmearingEnum.FERMI.value)
        except AttributeError:
            pass
        try:
            if self._parameters.smearing.mp:
                self._set_simple('ismear', IntSmearingEnum.MP.value * abs(int(self._parameters.smearing.mp)))
        except AttributeError:
            pass
        try:
            if self._parameters.smearing.tetra:
                self._set_simple('ismear', IntSmearingEnum.TETRA.value)
        except AttributeError:
            pass

    def set_icharg(self):  # noqa: MC0001
        """
        Set the flag to start from input charge density and keep it constant.

        See: https://www.vasp.at/wiki/index.php/ICHARG
        """
        try:
            if self._parameters.charge.from_wave:
                self._set_simple('icharg', ChargeEnum.WAVE.value)
        except AttributeError:
            pass
        try:
            if self._parameters.charge.from_charge:
                self._set_simple('icharg', ChargeEnum.CHARGE.value)
        except AttributeError:
            pass
        try:
            if self._parameters.charge.from_atomic:
                self._set_simple('icharg', ChargeEnum.ATOMIC.value)
        except AttributeError:
            pass
        try:
            if self._parameters.charge.from_potential:
                self._set_simple('icharg', ChargeEnum.POTENTIAL.value)
        except AttributeError:
            pass
        try:
            if self._parameters.charge.constant_charge:
                self._set_simple('icharg', ChargeEnum.CONSTANT_CHARGE.value)
        except AttributeError:
            pass
        try:
            if self._parameters.charge.constant_atomic:
                self._set_simple('icharg', ChargeEnum.CONSTANT_ATOMIC.value)
        except AttributeError:
            pass

    def set_lorbit(self):  # noqa: MC0001
        """
        Set the flag that controls the projectors/decomposition onto orbitals.

        See: https://www.vasp.at/wiki/index.php/LORBIT
        """
        self._set_wigner_seitz_radius()
        try:
            if self._parameters.bands.decompose_bands:
                if self._parameters.bands.decompose_wave:
                    # Issue a warning that one can only use either or
                    raise ValueError('Only projections/decompositions on the bands or the wave function are allowed.')
                wigner_seitz_radius = False
                try:
                    if abs(self._massage.rwigs[0]) > 1E-8:
                        wigner_seitz_radius = True
                except AttributeError:
                    pass
                if self._parameters.bands.decompose_auto:
                    self._set_simple('lorbit', OrbitEnum.NO_RWIGS_ATOM_LM_PHASE_AUTO.value)
                else:
                    try:
                        lm = self._parameters.bands.lm  # pylint: disable=invalid-name
                    except AttributeError:
                        lm = False  # pylint: disable=invalid-name
                    try:
                        phase = self._parameters.bands.phase
                    except AttributeError:
                        phase = False
                    lorbit = OrbitEnum.get_lorbit_from_combination(lm=lm, phase=phase, wigner_seitz_radius=wigner_seitz_radius).value
                    self._set_simple('lorbit', lorbit)
            else:
                try:
                    if self._parameters.bands.decompose_wave:
                        self._set_simple('lorbit', OrbitEnum.ATOM_LM_WAVE.value)
                except AttributeError:
                    pass
        except AttributeError:
            try:
                if self._parameters.bands.decompose_wave:
                    self._set_simple('lorbit', OrbitEnum.ATOM_LM_WAVE.value)
            except AttributeError:
                pass

    def _set_wigner_seitz_radius(self):
        """
        Set the Wigner Seitz radius that is used to project/decompose.

        See: https://www.vasp.at/wiki/index.php/RWIGS
        """
        try:
            wigner_seitz_radius = self._parameters.bands.wigner_seitz_radius
            # Check that it is defined as a list
            if isinstance(wigner_seitz_radius, list):
                if wigner_seitz_radius[0]:
                    self._set_simple('rwigs', wigner_seitz_radius)
            else:
                raise ValueError('The parameter wigner_seitz_radius should be supplied as a list of floats bigger than zero.')
        except AttributeError:
            pass

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


def inherit_and_merge_parameters(inputs):
    """
    Goes trough the inputs namespaces and the namespaces in the inputs.parameters and merge them.

    Note that parameters specified in the inputs.parameters will override what is supplied as workchain input,
    in case there is overlap.
    """
    parameters = AttributeDict()
    namespaces = ['electronic', 'bands', 'smearing', 'charge', 'relax', 'converge']
    for namespace in namespaces:  # pylint: disable=too-many-nested-blocks
        parameters[namespace] = AttributeDict()
        try:
            for key, item in inputs[namespace].items():
                if isinstance(item, DataFactory('array')):
                    # Only allow one array per input
                    if len(item.get_arraynames()) > 1:
                        raise IndexError(
                            'The input array with a key {} contains more than one array. Please make sure an input only contains one array.'
                            .format(key))
                    for array in item.get_arraynames():
                        parameters[namespace][key] = item.get_array(array)
                elif isinstance(item, DataFactory('dict')):
                    parameters[namespace][key] = item.get_dict()
                elif isinstance(item, DataFactory('list')):
                    parameters[namespace][key] = item.get_list()
                else:
                    parameters[namespace][key] = item.value
        except KeyError:
            pass

    # Now get the input parameters and update the dictionary. This means,
    # any supplied namespace in the parameters (i.e. inputs.parameters.somekey) will override what is supplied to the workchain
    # input namespace (i.e. inputs.somekey).
    try:
        # inputs might not have parameters, or parameters might be empty
        input_parameters = AttributeDict(inputs.parameters.get_dict())
    except AttributeError:
        input_parameters = {}

    # Now check that no loose keys are residing on the root of input_parameters, everything should be in
    # the vasp or aiida namespace
    #valid_keys = ['vasp', 'aiida']
    #if not list(input_parameters.keys()).sort() == valid_keys.sort():
    #    raise ValueError('Unsupported keys detected on parameter root. '
    #                     'Please make sure all keys reside inside the vasp or aiida namespace.')

    # We cannot use regular update here, as we only want to replace each key if it exists, if a key
    # contains a new dict we need to traverse that, hence we have a function to perform this update
    update_nested_dict(parameters, input_parameters)

    return parameters
