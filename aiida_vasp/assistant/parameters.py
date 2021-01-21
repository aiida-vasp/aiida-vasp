"""
Parameter related utils

-----------------------
Contains utils and definitions that are used together with the parameters.
"""
# pylint: disable=too-many-branches

import enum

from aiida.common.extendeddicts import AttributeDict
from aiida.plugins import DataFactory
from aiida_vasp.utils.extended_dicts import update_nested_dict

_BASE_NAMESPACES = ['electronic', 'smearing', 'charge', 'dynamics', 'bands', 'relax', 'converge']
_DEFAULT_OVERRIDE_NAMESPACE = 'incar'


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


class ParametersMassage():  # pylint: disable=too-many-instance-attributes
    """
    A class that contains all relevant massaging of the input parameters for VASP.

    The idea is that this class accepts the set input parameters from AiiDA (non code specifics), checks if any code specific
    parameters supplied are valid VASP input parameters (only rudimentary at this point, should also cross check and check types)
    and convert the AiiDA input parameters to VASP specific parameters. A set function needs to be developed for each VASP INCAR
    parameter that we want to set based on supplied AiiDA/AiiDA-VASP specific parameters. These set functions takes these parameters
    and converts it to VASP INCAR compatible tags. The parameter property should return ready to go parameters containing the
    default override namespace, the namespaces set in the `_set_extra_parameters` function and any additional namespaces
    that might have been set using the `additional_override_namespaces` setting entry in settings
    that can be supplied to the VaspWorkChain.

    The default override namespace (see `_DEFAULT_OVERRIDE_NAMESPACE`) should always be present when using this VASP plugin.
    If using additional plugins, one can for instance supply additional namespace override that can be used,
    depending on what is needed in those plugins and how you construct your workchains.
    """

    def __init__(self, parameters, unsupported_parameters=None, settings=None):
        self.exit_code = None

        # Check type of parameters and set
        self._parameters = check_inputs(parameters)

        # Check type of supplied unsupported_parameters and set
        self._unsupported_parameters = check_inputs(unsupported_parameters)

        # Check type of settings and set
        self._settings = check_inputs(settings)

        # Check setting for any possible supplied additional override namespace and set
        self._additional_override_namespaces = self._fetch_additional_override_namespaces()

        self._massage = AttributeDict()
        # Initialize any allowed override namespaces
        self._massage[_DEFAULT_OVERRIDE_NAMESPACE] = AttributeDict()
        for item in self._additional_override_namespaces:
            self._massage[item] = AttributeDict()

        self._check_valid_namespaces()
        # Load the valid INCAR parameters that are supported by VASP
        self._load_valid_params()
        # Establish set functions which will convert AiiDA/AiiDA-VASP specific parameters to VASP INCAR tags
        self._functions = ParameterSetFunctions(self._parameters, self._massage[_DEFAULT_OVERRIDE_NAMESPACE])
        # Convert general parameters to VASP specific ones, using the set functions
        self._set_vasp_parameters()
        # Override or set parameters that are supplied in the default override namespace (should be valid VASP INCAR tags)
        self._set_override_vasp_parameters()
        # Set any extra parameters not related to INCAR
        self._set_extra_vasp_parameters()
        # Set any additional override namespace that should just be forwarded
        self._set_additional_override_parameters()
        # Finally, we validate the INCAR parameters in order to prepare it for dispatch
        self._validate_vasp_parameters()

    def _check_valid_namespaces(self):
        """Check that we do not have namespaces on the input parameters that is unsupported."""
        for key in self._parameters.keys():
            if key not in list(_BASE_NAMESPACES + self._additional_override_namespaces + [_DEFAULT_OVERRIDE_NAMESPACE]):
                raise ValueError(f'The supplied namespace: {key} is not supported.')

    def _load_valid_params(self):
        """Import a list of valid parameters for VASP. This is generated from the manual."""
        from os import path  # pylint: disable=import-outside-toplevel
        from yaml import safe_load  # pylint: disable=import-outside-toplevel
        with open(path.join(path.dirname(path.realpath(__file__)), 'parameters.yml'), 'r') as file_handler:
            tags_data = safe_load(file_handler)
        self._valid_parameters = list(tags_data.keys())
        # Now add any unsupported parameter to the list
        for key, _ in self._unsupported_parameters.items():
            key = key.lower()
            if key not in self._valid_parameters:
                self._valid_parameters.append(key)

    def _fetch_additional_override_namespaces(self):
        """Check the settings for any additional supplied override namespace and return it."""
        try:
            override_namespaces = self._settings.additional_override_namespaces
        except AttributeError:
            override_namespaces = []

        return override_namespaces

    def _set_vasp_parameters(self):
        """Iterate over the valid parameters and call the set function associated with that parameter."""
        for key in self._valid_parameters:
            self._set(key)

    def _set_override_vasp_parameters(self):
        """Set the any supplied override parameters."""
        try:
            if self._parameters[_DEFAULT_OVERRIDE_NAMESPACE]:
                for key, item in self._parameters[_DEFAULT_OVERRIDE_NAMESPACE].items():
                    # Sweep the override input parameters (only care about the ones in the default override namespace)
                    # to check if they are valid VASP tags.
                    key = key.lower()
                    if self._valid_vasp_parameter(key):
                        # Add or override in the default override namespace
                        self._massage[_DEFAULT_OVERRIDE_NAMESPACE][key] = item
                    else:
                        break
        except KeyError:
            # The default override namespace might not be supplied (no override)
            pass

    def _set_extra_vasp_parameters(self):
        """
        Find if there are any extra parameters that are not part of the INCAR that needs to be set.

        One example is the dynamic namespace which handles for instance flags for selective dynamics.
        These flags are more connected to a calculation than a StructureData and thus it was necessary
        to make sure it was valid input to the VASP workchain.

        """
        try:
            if self._parameters.dynamics:
                self._massage.dynamics = AttributeDict()
            for key, item in self._parameters.dynamics.items():
                key = key.lower()
                if key in ['positions_dof']:
                    self._massage.dynamics[key] = item
                else:
                    break
        except AttributeError:
            pass

    def _set_additional_override_parameters(self):
        """Set any customized parameter namespace, including its content on the massaged container."""
        parameters_keys = self._parameters.keys()
        for item in self._additional_override_namespaces:
            if item in parameters_keys:
                # Only add if namespace exists in parameters
                self._massage[item] = AttributeDict(self._parameters[item])

    def _valid_vasp_parameter(self, key):
        """Make sure a key are recognized as a valid VASP input parameter."""
        if key not in self._valid_parameters:
            raise ValueError(f'The supplied key: {key} is not a support VASP parameter.')

        return True

    def _validate_vasp_parameters(self):
        """Make sure all the massaged values are recognized as valid VASP input parameters."""
        for key in self._massage[_DEFAULT_OVERRIDE_NAMESPACE]:
            key = key.lower()
            if not self._valid_vasp_parameter(key):
                break

    def _set(self, key):
        """Call the necessary function to set each parameter."""
        try:
            getattr(self._functions, 'set_' + key)()
        except AttributeError:
            # We have no setter function for the valid key, meaning there is no general parameter that is linked
            # to this key (INCAR tag). These tags have to be supplied in the default override namespace for now.
            pass

    @property
    def parameters(self):
        """Return the massaged parameter set ready to go in VASP format."""
        return self._massage


class ParameterSetFunctions():
    """Container for the set functions that converts an AiiDA parameters to a default override specific one."""

    def __init__(self, parameters, incar):
        self._parameters = parameters
        self._incar = incar

    def set_encut(self):
        """
        Set which plane wave cutoff to use.

        See https://www.vasp.at/wiki/index.php/ENCUT
        """

        try:
            self._incar.encut = self._parameters.electronic.pwcutoff
        except AttributeError:
            pass

    def set_ibrion(self):
        """
        Set which algorithm to use for ionic movements.

        See: https://www.vasp.at/wiki/index.php/IBRION
        """

        if self._relax():
            try:
                if self._parameters.relax.algo == 'cg':
                    self._incar.ibrion = RelaxAlgoEnum.IONIC_RELAXATION_CG.value
                elif self._parameters.relax.algo == 'rd':
                    self._incar.ibrion = RelaxAlgoEnum.IONIC_RELAXATION_RMM_DIIS.value
                else:
                    raise ValueError(f'The supplied relax.algo: {self._parameters.relax.algo} is not supported')
            except AttributeError:
                pass

    def set_ediffg(self):
        """
        Set the cutoff to use for relaxation.

        See: https://www.vasp.at/wiki/index.php/EDIFFG
        """

        if not self._relax():
            # This flag is only valid if you have enabled relaxation
            return
        energy_cutoff = False
        try:
            self._incar.ediffg = self._parameters.relax.energy_cutoff
            energy_cutoff = True
        except AttributeError:
            pass
        try:
            self._incar.ediffg = -abs(self._parameters.relax.force_cutoff)
            if energy_cutoff:
                raise ValueError('User supplied both a force and an energy cutoff for the relaxation. Please select.')
        except AttributeError:
            pass

    def set_nsw(self):
        """
        Set the number of ionic steps to perform.

        See: https://www.vasp.at/wiki/index.php/NSW
        """

        if self._relax():
            try:
                self._set_simple('nsw', self._parameters.relax.steps)
            except AttributeError:
                pass

    def set_isif(self):
        """
        Set relaxation mode according to the chosen degrees of freedom.

        See: https://www.vasp.at/wiki/index.php/ISIF
        """

        if self._relax():
            positions = self._parameters.get('relax', {}).get('positions', False)
            shape = self._parameters.get('relax', {}).get('shape', False)
            volume = self._parameters.get('relax', {}).get('volume', False)
            try:
                self._incar.isif = RelaxModeEnum.get_isif_from_dof(positions=positions, shape=shape, volume=volume).value
            except AttributeError:
                pass

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
                    if abs(self._incar.rwigs[0]) > 1E-8:
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
        return bool(self._parameters.get('relax', {}).get('positions') or \
            self._parameters.get('relax', {}).get('shape') or \
            self._parameters.get('relax', {}).get('volume'))

    def _set_simple(self, target, value):
        """Set basic parameter."""
        try:
            self._incar[target] = value
        except AttributeError:
            pass


def check_inputs(supplied_inputs):
    """Check that the inputs are of some correct type and returned as AttributeDict."""
    inputs = None
    if supplied_inputs is None:
        inputs = AttributeDict()
    else:
        if isinstance(supplied_inputs, DataFactory('dict')):
            inputs = AttributeDict(supplied_inputs.get_dict())
        elif isinstance(supplied_inputs, dict):
            inputs = AttributeDict(supplied_inputs)
        elif isinstance(supplied_inputs, AttributeDict):
            inputs = supplied_inputs
        else:
            raise ValueError(f'The supplied type {type(inputs)} of inputs is not supported. Supply a dict, Dict or an AttributeDict.')

    return inputs


def inherit_and_merge_parameters(inputs):
    """
    Goes trough the inputs namespaces and the namespaces in the inputs.parameters and merge them.

    Note that parameters specified in the inputs.parameters will override what is supplied as workchain input,
    in case there is overlap.
    """
    parameters = AttributeDict()
    namespaces = _BASE_NAMESPACES

    # We start with a clean parameters and first set the allowed namespaces and its content from the inputs of the workchain
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

    # Then obtain the inputs.parameters.
    # Here we do not do any checks for valid parameters, that is done later when reaching the ParameterMassager.
    try:
        input_parameters = AttributeDict(inputs.parameters.get_dict())
    except AttributeError:
        # Inputs might not have parameters
        input_parameters = AttributeDict()

    # Now the namespace and content of the workchain inputs and the inputs.parameters are merged.
    # Any supplied namespace in the parameters (i.e. inputs.parameters.somekey) will override what
    # is supplied to the workchain input namespace (i.e. inputs.somekey).
    # We cannot use regular update here, as we only want to replace each key if it exists, if a key
    # contains a new dict we need to traverse that, hence we have a function to perform this update.
    update_nested_dict(parameters, input_parameters)

    return parameters
