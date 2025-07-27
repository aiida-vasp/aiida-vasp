"""
Parser module for composing output of a VASP calculation.

The simplified parser outputs the following nodes:

1. A `misc` node that stores simple summary information such as total energies,
    total run times, any warnings issues during the calculation and if the calculation
    was finished.
2. A `arrays` node that for storing miscellaneous quantities that are arrays by nature and typically
    have a large size.
3. A `trajectory` node for storing the trajectory of geometry optimisation and AIMD.
4. A `bands` node for storing the band structure.
5. A `dos` node for storing the density of states.
6. Other nodes for storing other relevant quantities such as the born effective charges.


Main difference from the previous version

1. `pydantic` is used to validate the parser settings *at submission* time.
2. When parsing retrieved data, we take a 'parse as much as possible' approach -
   the quantities are always parsed if available, and only excluded during the stage of composing the output nodes.
   This is simpler from the previous 'parse only when needed' approach, where multiple checks have be done to work
   out which parser to call and which quantities to include.
3. All parser logic is contained in a single class and can be extended by updating content parsers and modify the
   default settings and update/add the `_compose_xx` methods.
"""

from typing import Any, Dict, List

import numpy as np
from aiida import orm
from aiida.engine import ExitCode
from aiida.parsers.parser import Parser
from pydantic import Field

from aiida_vasp.data.chargedensity import ChargedensityData
from aiida_vasp.data.wavefun import WavefunData
from aiida_vasp.parsers.content_parsers import *
from aiida_vasp.utils.opthold import OptionContainer


class ParserError(RuntimeError):
    pass


class QuantityMissingError(ParserError):
    pass


class RequiredQuantityMissingError(ParserError):
    pass


class MissingFileError(ParserError):
    pass


DEFAULT_EXCLUDED_QUANTITIES = (
    'elastic_moduli',
    'symmetries',
    'parameters',  # The parameters used for the calculation
)

DEFAULT_EXCLUDED_NODE = ('bands', 'dos', 'kpoints', 'trajectory', 'energies', 'wavecar', 'chgcar', 'projectors')

DEFAULT_REQUIRED_QUANTITIES = ('run_status', 'run_stats')

DEFAULT_FILE_MAPPING = {
    'vasprun.xml': 'vasprun.xml',
    'vasp_output': 'vasp_output',
    'OUTCAR': 'OUTCAR',
    'CONTCAR': 'CONTCAR',
    'CHGCAR': 'CHGCAR',
    'IBZKPT': 'IBZKPT',
}
MISC_QUANTITIES = (
    'total_energies',
    'notifications',
    'run_status',
    'run_stats',
    'version',
    'forces',
    'stress',
    'site_magnetization',
    'band_properties',
    'elastic_moduli',
    'symmetries',
    'fermi_level',
    'band_properties',
    'magnetization',
    'parameters',
)

ALLOW_EMPTY = ('notifications',)

# Quantities that should be stored inside separate nodes
STANDALONE_ARRAY_QUANTITIES = {
    'born_charges': 'vasprun.xml',
    'dielectrics': 'vasprun.xml',
    'dynmat': 'vasprun.xml',
    'hessian': 'vasprun.xml',
    'projectors': 'vasprun.xml',
    'energies': 'vasprun.xml',
}


class ParserSettingsConfig(OptionContainer):
    """
    Settings for the VASP parser.
    """

    include_quantity: List[str] = Field(description='Properties to include', default_factory=lambda: [])
    exclude_quantity: List[str] = Field(description='Quantities to be excluded', default_factory=lambda: [])
    required_quantity: List[str] = Field(
        description='Quantities that most be present', default_factory=lambda: list(DEFAULT_REQUIRED_QUANTITIES)
    )
    include_node: List[str] = Field(description='Output node to include', default_factory=lambda: [])
    exclude_node: List[str] = Field(description='Output node to exclude', default_factory=lambda: [])
    file_mapping: Dict[str, str] = Field(
        description='Mapping of file names to quantities', default_factory=lambda: dict(DEFAULT_FILE_MAPPING)
    )
    kpoints_from_ibzkpt: bool = False
    check_completeness: bool = True
    electronic_step_energies: bool = False
    energy_type: List[str] = Field(
        description='Energy types to include', default_factory=lambda: ['energy_extrapolated']
    )
    keep_stream_history: bool = Field(
        description='Whether to keep the history of all notifications in the parsed stream (stdout)', default=False
    )
    ignore_notification_errors: bool = Field(
        description='Whether to ignore errors in the notifications parsed from vasp_output', default=False
    )
    critical_notification_errors: List[str] = Field(
        description='Critical stream errors to raise',
        default_factory=lambda: [
            'brmix',
            'edwave',
            'brmix',
            'cnormn',
            'denmp',
            'dentet',
            'edddav_zhegv',
            'eddrmm_zhegv',
            'edwav',
            'fexcp',
            'fock_acc',
            'non_collinear',
            'not_hermitian',
            'pzstein',
            'real_optlay',
            'rhosyg',
            'rspher',
            'set_indpw_full',
            'sgrcon',
            'no_potimm',
            'magmom',
            'bandocc',
        ],
    )
    critical_objects: List[str] = Field(
        description='Critical objects to be present', default_factory=lambda: ['vasprun.xml', 'OUTCAR']
    )
    check_errors: bool = Field(description='Whether to check for errors in calculation', default=True)
    check_ionic_convergence: bool = Field(
        description='Whether to check for convergence during the relaxation based on the INCAR settings', default=True
    )
    omit_structure: bool = Field(
        description='Whether to omit the structure node from the output if no ionic movement', default=True
    )


class VaspParser(Parser):
    """Class for parsing VASP output files and storing the results in AiiDA."""

    def __init__(self, node: orm.CalcJobNode) -> None:
        """
        Initialize the Parser instance
        """
        super(VaspParser, self).__init__(node)
        # Create the containers
        self.user_config = None
        self.quantities_each: Dict[str, Any] = {}
        self.errored_quantities: Dict[str, Any] = {}
        self.errored_parsers: Dict[str, Any] = {}
        self.parser_notifications: Dict[str, List[str]] = {}
        self.retrieve_object_names: List[str] = []
        self.quantities_to_exclude: List[str] = []
        self.nodes_to_exclude: List[str] = []

    def _init_user_settings(self) -> ParserSettingsConfig:
        """Initialize the settings from the inputs."""
        if 'settings' in self.node.inputs:
            user_config: ParserSettingsConfig = ParserSettingsConfig(
                **self.node.inputs.settings.get_dict().get('parser_settings', {})
            )
        else:
            user_config = ParserSettingsConfig()
        # Initialize the containers
        self.user_config = user_config
        return user_config

    def _get_quantities_to_parse(self) -> ExitCode | None:
        """Return the list of quantities to parse."""
        # Apply the modifiers
        user_config = self.user_config
        quantities_to_exclude = [key for key in DEFAULT_EXCLUDED_QUANTITIES if key not in user_config.include_quantity]
        quantities_to_exclude += user_config.exclude_quantity
        nodes_to_exclude = [key for key in DEFAULT_EXCLUDED_NODE if key not in user_config.include_node]
        nodes_to_exclude += user_config.exclude_node
        self.quantities_to_exclude = quantities_to_exclude
        self.nodes_to_exclude = nodes_to_exclude

        # Check for critical missing objects
        self.retrieve_object_names = self.retrieved.list_object_names()
        missing = False
        for name, _ in user_config.file_mapping.items():
            if name in user_config.critical_objects and name not in self.retrieve_object_names:
                missing = True
        if missing is True:
            return self.exit_codes.ERROR_CRITICAL_MISSING_OBJECT

    def _post_process_quantities(self) -> ExitCode | None:
        """Post-process the parsed quantities."""

        # Warn about errored/missing quantities and parsers

        if self.errored_quantities:
            self.logger.warning(
                f'The following quantities cannot be parsed due to errors: {", ".join(self.errored_quantities)}'
            )
        if self.errored_parsers:
            self.logger.warning(
                f'The following parsers cannot be instantiated due to: {", ".join(self.errored_parsers)}'
            )

        # Remove the quantities
        for name, parsed_quantities in self.quantities_each.items():
            for sub_key in list(parsed_quantities.keys()):
                if sub_key in self.quantities_to_exclude:
                    del parsed_quantities[sub_key]

        # Check in required quantities are present
        missing_required = []
        for name in self.user_config.required_quantity:
            exists = False
            for _, value in self.quantities_each.items():
                if value.get(name) is not None:
                    exists = True
                    break
            if exists is False:
                missing_required.append(name)
        if missing_required:
            return self.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY.format(quantity=','.join(missing_required))

    def parse(self, **kwargs: Any) -> ExitCode | None:
        """
        Parse outputs, store results in database.
        """
        user_config = self._init_user_settings()

        exit_code = self._get_quantities_to_parse()
        if exit_code is not None:
            return exit_code

        # Parse the files
        def parse_and_add(
            name: str,
            parser_cls: Any,
            required: bool = True,
            open_mode: str = 'r',
            content_parser_settings: dict | None = None,
        ) -> None:
            """Parse the target file and add the result to the quantities_each dictionary"""
            resolved_name = user_config.file_mapping[name]
            if resolved_name in self.retrieve_object_names:
                with self.retrieved.open(resolved_name, open_mode) as handler:
                    try:
                        parser: BaseFileParser = parser_cls(handler=handler, settings=content_parser_settings)
                    except Exception as error:
                        self.errored_parsers[name] = error
                        return
                    if parser.parser_notifications:
                        self.parser_notifications.update(parser.parser_notifications)
                    self.quantities_each[name], errored = parser.get_all_quantities()
                    self.errored_quantities.update(errored)
            elif user_config.check_completeness is True and required is True:
                raise MissingFileError(f'{resolved_name} is missing in the retrieved folder.')

        parse_and_add(
            'vasprun.xml',
            VasprunParser,
            required=True,
            open_mode='rb',
            content_parser_settings={
                'electronic_step_energies': user_config.electronic_step_energies,
                'energy_type': user_config.energy_type,
                'stream_history': user_config.keep_stream_history,
            },
        )
        parse_and_add('OUTCAR', OutcarParser, required=True)
        parse_and_add('vasp_output', StreamParser, required=True)
        parse_and_add('CONTCAR', PoscarParser, required=True)

        if user_config.kpoints_from_ibzkpt:
            parse_and_add('IBZKPT', KpointsParser, required=True)

        exit_code = self._post_process_quantities()
        if exit_code is not None:
            return exit_code

        return self._create_outputs()

    def _create_outputs(self) -> ExitCode | None:
        """Create the output nodes"""
        # Create the outputs
        self._failed_to_compose = {}

        # Call the _compose_xx methods to create the output nodes
        for method_name in [item for item in self.__dir__() if item.startswith('_compose_')]:
            name = method_name.replace('_compose_', '')
            if name in self.nodes_to_exclude:
                continue
            node_or_dict = None
            try:
                node_or_dict = getattr(self, '_compose_' + name)(self.quantities_each)
            except (QuantityMissingError, KeyError, ValueError, TypeError, AttributeError) as error:
                self._failed_to_compose[name] = error
                self.logger.warning(f'Failed to compose {name} node: {error}')
                continue
            if isinstance(node_or_dict, orm.Data):
                self.out(name, node_or_dict)
            elif isinstance(node_or_dict, dict):
                for key, value in node_or_dict.items():
                    self.out(key, value)
        if (
            any(name in self.user_config.include_node for name in self._failed_to_compose)
            and self.user_config.check_completeness is True
        ):
            return self.exit_codes.ERROR_NOT_ABLE_TO_CREATE_NODE.format(nodes=', '.join(self._failed_to_compose.keys()))
        # Check for errors
        if self.user_config.check_errors is True:
            error = self._check_vasp_errors(self.parser_notifications)
            return error

    def _compose_misc(self, quantities_each: dict[str, Any]) -> orm.Dict:
        """Compose the `misc` output node"""

        out_dict = {}
        gather_quantities(quantities_each, self.user_config.file_mapping['vasprun.xml'], out_dict, MISC_QUANTITIES)
        gather_quantities(quantities_each, self.user_config.file_mapping['OUTCAR'], out_dict, MISC_QUANTITIES)
        gather_quantities(quantities_each, self.user_config.file_mapping['vasp_output'], out_dict, MISC_QUANTITIES)
        # Filter field with all empty container
        out_dict = {key: value for key, value in out_dict.items() if not is_all_empty(value) or key in ALLOW_EMPTY}
        return orm.Dict(dict=out_dict)

    def _compose_structure(self, quantities_each: dict[str, Any]) -> orm.StructureData | None:
        """Compose the `structure` output node"""

        data = None
        # Omit output structure if not doing ionic relaxation
        # Better to inspect the parameters recorded directly inside the vasprun.xml
        if 'parameters' in self.node.inputs:
            incar_dict = {key.lower(): value for key, value in self.node.inputs.parameters.get_dict().items()}
            if (
                incar_dict.get('ibrion', -1) < 0 or incar_dict.get('nsw', 0) <= 0
            ) and self.user_config.omit_structure is True:
                self.logger.info('No ionic movement detected, omitting the structure output node.')
                return None

        if 'vasprun.xml' in quantities_each:
            data = quantities_each['vasprun.xml'].get('structure')
        if data is None:
            data = quantities_each.get('CONTCAR', {}).get('structure')
        if data is None:
            raise QuantityMissingError()
        return get_structure_node(data)

    def _compose_wavecar(self, quantities_each: dict[str, Any]) -> None:
        """Compose the `wavecar` output node"""

        # Check if WAVECAR is present in the retrieved folder
        if 'WAVECAR' in self.retrieve_object_names:
            with self.retrieved.base.repository.open('WAVECAR', 'rb') as handler:
                self.outputs['wavecar'] = WavefunData(file=handler, filename='WAVECAR')
        else:
            self.logger.warning('WAVECAR is not present in the retrieved folder.')

    def _compose_chgcar(self, quantities_each: dict[str, Any]) -> None:
        """Compose the `chgcar` output node"""

        # Check if WAVECAR is present in the retrieved folder
        if 'CHGCAR' in self.retrieve_object_names:
            with self.retrieved.base.repository.open('CHGCAR', 'rb') as handler:
                self.outputs['chgcar'] = ChargedensityData(file=handler, filename='CHGCAR')
        else:
            self.logger.warning('CHGCAR is not present in the retrieved folder.')

    def _compose_arrays(self, quantities_each: dict[str, Any]) -> dict[str, orm.ArrayData]:
        """Generate the generic `arrays` output node"""
        out_arrays = {}

        # Compose the standalone arrays - each corresponds to a single quantity
        for name, file_name in STANDALONE_ARRAY_QUANTITIES.items():
            if name in self.nodes_to_exclude:
                continue
            array_node = self._make_standalone_array(quantities_each, name, file_name)
            if array_node is not None:
                out_arrays[name] = array_node
        return out_arrays

    def _make_standalone_array(
        self, quantities_each: dict[str, Any], name: str, file_name: str = 'vasprun.xml'
    ) -> orm.ArrayData | None:
        """Compose the `dielectrics` output node"""
        # The output can be an array or a dictionary of arrays - both cases should be handled
        arrays_or_dict = quantities_each.get(file_name, {}).get(name)
        # Avoid creating empty arrays
        if isinstance(arrays_or_dict, dict) and len(arrays_or_dict) > 0:
            arrays_or_dict = {key: value for key, value in arrays_or_dict.items() if value is not None}
            return orm.ArrayData(arrays_or_dict)
        elif isinstance(arrays_or_dict, (np.ndarray, list)) and len(arrays_or_dict) > 0:
            return orm.ArrayData({name: arrays_or_dict})
        return None

    def _compose_kpoints(self, quantities_each: dict[str, Any]) -> orm.KpointsData:
        """Compose the `kpoints` output node"""
        kpoints_data = None
        if self.user_config.kpoints_from_ibzkpt is True:
            kpoints_data = quantities_each['IBZKPT']['kpoints']
        elif 'vasprun.xml' in quantities_each:
            kpoints_data = quantities_each['vasprun.xml'].get('kpoints')

        if kpoints_data is not None:
            node = orm.KpointsData()
            if kpoints_data['mode'] == 'explicit':
                node.set_kpoints(
                    kpoints_data['points'], weights=kpoints_data['weights'], cartesian=kpoints_data['cartesian']
                )
            elif kpoints_data['mode'] == 'automatic':
                node.set_kpoints_mesh(kpoints_data['divisions'], offset=kpoints_data['shifts'])
            else:
                raise ValueError(f'Unknown kpoints mode {kpoints_data["mode"]}')
            # Record the cell for which the kpoints are defined for
            node.set_cell(quantities_each['vasprun.xml']['structure']['unitcell'])
            return node
        raise QuantityMissingError('No valid kpoints data to use')

    def _compose_trajectory(self, quantities_each: dict[str, Any]) -> orm.TrajectoryData | None:
        """Compose the `trajectory` output"""

        if 'vasprun.xml' in quantities_each:
            node = orm.TrajectoryData()
            traj_data = quantities_each['vasprun.xml'].get('trajectory')
            # No need to carry on if there are no trajectory data
            if traj_data is None or len(traj_data) == 0:
                return None
            for key, value in traj_data.items():
                if key == 'symbols':
                    node.base.attributes.set(key, value)
                elif value.dtype.hasobject:
                    self.logger.warning(f'Cannot set array {key}: {value} in TrajectoryData as it is not numerical.')
                else:
                    node.set_array(key, value)
            for key, value in quantities_each['vasprun.xml']['energies'].items():
                node.set_array(key, value)
            return node
        return None

    def _compose_bands(self, quantities_each: dict[str, Any]) -> orm.BandsData:
        """Compose the `band` node"""
        if 'vasprun.xml' in quantities_each:
            deigen = quantities_each['vasprun.xml']['eigenvalues']
            docc = quantities_each['vasprun.xml']['occupancies']
            if 'total' in deigen:
                eigenvalues = np.array(deigen['total'])
                occupancies = np.array(docc['total'])
            else:
                eigenvalues = np.array([deigen['up'], deigen['down']])
                occupancies = np.array([docc['up'], docc['down']])
            node = orm.BandsData()
            kpoints = self._compose_kpoints(quantities_each)
            node.set_kpointsdata(kpoints)
            node.set_bands(eigenvalues, occupations=occupancies)

            # Record the Fermi level if available
            node.base.attributes.set('fermi_level', quantities_each['vasprun.xml'].get('fermi_level'))
            node.base.attributes.set('efermi', quantities_each['vasprun.xml'].get('fermi_level'))
            node.set_cell(quantities_each['vasprun.xml']['structure']['unitcell'])
            return node

    def _compose_dos(self, quantities_each: dict[str, Any]) -> orm.ArrayData | None:
        """Compose the `dos` node"""
        arrays_dict = {}
        if 'vasprun.xml' in quantities_each:
            gather_quantities(quantities_each, 'dos', arrays_dict, ['dos'], flatten_dict=True)
        if arrays_dict:
            node = orm.ArrayData(arrays_dict['dos'])
            return node

    def _check_vasp_errors(self, parser_notifications: dict[str, Any]) -> ExitCode | None:
        """
        Detect simple vasp execution problems and returns the exit_codes to be set
        """
        quantities = {}
        for key, value in self.quantities_each.items():
            for key_, value_ in value.items():
                quantities[key_] = value_
        if 'run_status' not in quantities:
            return self.exit_codes.ERROR_DIAGNOSIS_OUTPUTS_MISSING
        run_status = quantities['run_status']

        try:
            # We have an overflow in the XML file which is critical, but not reported by VASP in
            # the standard output, so checking this here.
            if parser_notifications.get('vasprun_xml_overflow'):
                return self.exit_codes.ERROR_OVERFLOW_IN_XML
        except AttributeError:
            pass

        # Return errors related to execution and convergence problems.
        # Note that the order is important here - if a calculation is not finished, we cannot
        # comment on wether properties are converged are not.
        if run_status['finished'] is False:
            return self.exit_codes.ERROR_DID_NOT_FINISH

        if run_status['electronic_converged'] is False:
            return self.exit_codes.ERROR_ELECTRONIC_NOT_CONVERGED

        # Check the ionic convergence issues
        if run_status['ionic_converged'] is False:
            if self.user_config.check_ionic_convergence is True:
                return self.exit_codes.ERROR_IONIC_NOT_CONVERGED
            self.logger.warning('The ionic relaxation is not converged, but the calculation is treated as successful.')

        # Check for the existence of critical warnings
        if 'notifications' in quantities:
            notifications = quantities['notifications']
            ignore_all = self.user_config.ignore_notification_errors
            if not ignore_all:
                composer = NotificationComposer(
                    notifications,
                    quantities['run_status'],
                    self.node.inputs,
                    self.exit_codes,
                    critical_notifications=self.user_config.critical_notification_errors,
                )
                exit_code = composer.compose()
                if exit_code is not None:
                    return exit_code
        else:
            self.logger.warning('WARNING: missing notification output for VASP warnings and errors.')

        return None


def gather_quantities(
    quantities_each: dict[str, Any], namespace: str, dst: dict[str, Any], fields: list[str], flatten_dict: bool = False
) -> None:
    """
    Gather quantities and put them into the target dictionary
    """
    for key, value in quantities_each.get(namespace, {}).items():
        if key in fields:
            if isinstance(value, dict) and flatten_dict:
                # flatten the dictionary - prepend the key with the name of the quantity
                for key2, value2 in value.items():
                    dst[key + '_' + key2] = value2
            else:
                dst[key] = value


class NotificationComposer:
    """Compose errors codes based on the notifications"""

    def __init__(
        self,
        notifications: list[dict[str, Any]],
        run_status: dict[str, Any],
        inputs: dict[str, Any],
        exit_codes: Any,
        critical_notifications: list[str],
    ) -> None:
        """
        Composed error codes based on the notifications

        Some of the errors need to have additional properties inspected before they can be emitted,
        as they might be trigged in a harmless way.

        To add new checkers, one needs to implement a property with the name of the error for this class and
        contains the code for checking. This property should return the exit_code or return None. The property
        is inspected if its name is in the list critical notifications.
        """
        self.notifications = notifications
        self.notifications_dict = {item['name']: item['message'] for item in self.notifications}
        self.run_status = run_status
        self.inputs = inputs
        self.exit_codes = exit_codes
        self.critical_notifications = critical_notifications

    def compose(self) -> ExitCode | None:
        """
        Compose the exit codes

        Returns None if no exit code should be emitted, otherwise emit the error code.
        """

        for critical in self.critical_notifications:
            # Check for any special handling
            if hasattr(self, critical):
                output = getattr(self, critical)
                if output:
                    return output
            # No special handling, just check if it exists
            elif critical in self.notifications_dict:
                return self.exit_codes.ERROR_VASP_CRITICAL_ERROR.format(error_message=self.notifications_dict[critical])
        return None

    @property
    def brmix(self) -> ExitCode | None:
        """Check if BRMIX should be emitted"""
        if 'brmix' not in self.notifications_dict:
            return None

        # If NELECT is set explicitly for the calculation then this is not an critical error
        if 'parameters' in self.inputs and 'nelect' in self.inputs['parameters'].get_dict():
            return None

        return self.exit_codes.ERROR_VASP_CRITICAL_ERROR.format(error_message=self.notifications_dict['brmix'])

    @property
    def edddav_zhegv(self) -> ExitCode | None:
        """Check if EDDDAV call to ZHEGV should be emitted. Sometimes it has converged."""
        if 'edddav_zhegv' not in self.notifications_dict:
            return None

        if self.run_status['electronic_converged']:
            return None

        return self.exit_codes.ERROR_VASP_CRITICAL_ERROR.format(error_message=self.notifications_dict['edddav_zhegv'])

    @property
    def eddrmm_zhegv(self) -> ExitCode | None:
        """Check if EDDRMM call to ZHEGV should be emitted. Sometimes it has converged."""
        if 'eddrmm_zhegv' not in self.notifications_dict:
            return None

        if self.run_status['electronic_converged']:
            return None

        return self.exit_codes.ERROR_VASP_CRITICAL_ERROR.format(error_message=self.notifications_dict['eddrmm_zhegv'])


def get_structure_node(structure_dict: dict[str, Any]) -> orm.StructureData:
    """Compose a structure node from the dictionary output by the parser"""
    node = orm.StructureData()
    node.set_cell(structure_dict['unitcell'])
    for site in structure_dict['sites']:
        node.append_atom(position=site['position'], symbols=site['symbol'], name=site['kind_name'])
    return node


def is_all_empty(obj: dict | list) -> bool:
    """Check if all elements of a dictionary or list are empty"""
    if isinstance(obj, dict):
        if len(obj) == 0:
            return True
        else:
            return all(is_all_empty(value) for value in obj.values())
    elif isinstance(obj, list):
        if len(obj) == 0:
            return True
        else:
            return all(is_all_empty(value) for value in obj)
    else:
        return False
