"""
Parser module for composing output of a VASP calculation.

The simplified parser outputs the following nodes:

1. A `misc` node that stores simply summary information such as total energies,
    total run times, any warnings issues during the calculation and if the calculation
    was finished.
2. A `arrays` node that for storing the properties that are arrays by nature
    such as the dielectric function, hessian matrix as well as the energies of each
    SCF cycle.
3. A `trajectory` node for storing the trajectory of geometry optimisation and AIMD.
4. A `band` node for storing the band structure.
5. A `dos` node for storing the density of states.

Main difference from the previous version
1. Individual properties are not stored in the separate node in order to reduce the number of nodes created. Scalar and
    small arrays are stored in misc (final energies and forces), and larger none-standard output arrays are stored in
    the `array` output node.
2. Standard outputs each as the bands and dos are stored in dedicated nodes.
3. `pydantic` is used to validate the parser settings *at submission* time.
4. When parsing retrieved data, we take a 'parse as much as possible' approach -
    the quantities are always parsed if available, and only excluded during the stage of composing the output nodes.
    This is simpler from the previous 'parse only when needed' approach, where multiple checks have be done to work
    out which parser to call and which quantities to include.
"""

from typing import Dict, List

import numpy as np
from aiida import orm
from aiida.parsers.parser import Parser
from pydantic import Field

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


DEFAULT_QUANTITIES = (
    'total_energies',
    'maximum_stress',
    'maximum_force',
    'notifications',
    'run_status',
    'run_stats',
    'version',
    'band_properties',
)

DEFAULT_EXCLUDED_QUANTITIES = (
    'energies',
    'trajectory',
    'kpoints',
    'chgcar',
    'wavecar',
    'projectors',
    'charge_density',
    'magnetization_density',
    'elastic_moduli',
    'symmetries',
    'band',
)

DEFAULT_EXCLUDED_NODE = tuple(['bands', 'dos', 'kpoints'])

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
    'maximum_stress',
    'maximum_force',
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
)


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
        description='Whether to keep the history of all notifications in the parsed stream' ' (stdout)', default=False
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


class VaspParser(Parser):
    """Class for parsing VASP output files and storing the results in AiiDA."""

    def __init__(self, node):
        """
        Initialize the Parser instance
        """
        super(VaspParser, self).__init__(node)

    def parse(self, **kwargs):
        """
        Parse outputs, store results in database.
        """

        if 'settings' in self.node.inputs:
            user_config: ParserSettingsConfig = ParserSettingsConfig(
                **self.node.inputs.settings.get_dict().get('parser_settings', {})
            )
        else:
            user_config = ParserSettingsConfig()
        self.user_config = user_config

        quantities_each = {}

        # Apply the modifiers
        quantities_to_exclude = [key for key in DEFAULT_EXCLUDED_QUANTITIES if key not in user_config.include_quantity]
        quantities_to_exclude += user_config.exclude_quantity
        nodes_to_exclude = [key for key in DEFAULT_EXCLUDED_NODE if key not in user_config.include_node]
        nodes_to_exclude += user_config.exclude_node

        # Check for critical missing objects
        retrieve_object_names = self.retrieved.list_object_names()
        missing = False
        for name, value in user_config.file_mapping.items():
            if name in user_config.critical_objects and name not in retrieve_object_names:
                missing = True
        if missing is True:
            return self.exit_codes.ERROR_CRITICAL_MISSING_OBJECT

        parser_notifications = {}
        errored_quantities = {}
        errored_parsers = {}

        self.errored_quantities = errored_quantities
        self.errored_parsers = errored_parsers

        # Parse the files
        def parse_and_add(name, parser_cls, required=True, open_mode='r', content_parser_settings=None):
            """Parse the target file and add the result to the quantities_each dictionary"""
            resolved_name = user_config.file_mapping[name]
            if resolved_name in retrieve_object_names:
                with self.retrieved.open(resolved_name, open_mode) as handler:
                    try:
                        parser: BaseFileParser = parser_cls(handler=handler, settings=content_parser_settings)
                    except Exception as error:
                        errored_parsers[name] = error
                        return
                    if parser.parser_notifications:
                        parser_notifications.update(parser.parser_notifications)
                    quantities_each[name], errored = parser.get_all_quantities()
                    errored_quantities.update(errored)
            elif user_config.check_completeness is True and required is True:
                raise MissingFileError(f'{resolved_name} is missing in the retrieved folder.')

        if errored_quantities:
            self.logger.warning(
                'The following quantities cannot be parsed due to errors:' f' {", ".join(errored_quantities)}'
            )
        if errored_parsers:
            self.logger.warning('The following parsers cannot be instantiated due to:' f' {", ".join(errored_parsers)}')

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

        if any(x not in quantities_to_exclude for x in ('charge_density', 'magnetization_density')):
            parse_and_add('CHGCAR', ChgcarParser, required=True)

        if user_config.kpoints_from_ibzkpt:
            parse_and_add('IBZKPT', ChgcarParser, required=True)
        self._quantities_each = quantities_each

        # Remove the quantities
        for name, parsed_quantities in quantities_each.items():
            for sub_key in list(parsed_quantities.keys()):
                if sub_key in quantities_to_exclude:
                    del parsed_quantities[sub_key]

        # Check in required quantities are present
        missing_required = []
        for name in user_config.required_quantity:
            exists = False
            for _, value in quantities_each.items():
                if value.get(name) is not None:
                    exists = True
                    break
            if exists is False:
                missing_required.append(name)
        if missing_required:
            return self.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY.format(quantity=','.join(missing_required))

        # Create the outputs
        self._failed_to_compose = {}
        for name in ['misc', 'structure', 'trajectory', 'kpoints', 'arrays', 'bands', 'dos']:
            if name in nodes_to_exclude:
                continue
            node = None
            try:
                node = getattr(self, '_compose_' + name)(quantities_each)
            except (QuantityMissingError, KeyError, ValueError, TypeError) as error:
                self._failed_to_compose[name] = error
                self.logger.warning(f'Failed to compose {name} node: {error}')
                continue
            if node is not None:
                self.out(name, node)
        if (
            any(name in user_config.include_node for name in self._failed_to_compose)
            and user_config.check_completeness is True
        ):
            return self.exit_codes.ERROR_NOT_ABLE_TO_CREATE_NODE.format(nodes=', '.join(self._failed_to_compose.keys()))
        # Check for errors
        if user_config.check_errors is True:
            error = self._check_vasp_errors(parser_notifications)
            return error

    def _compose_misc(self, quantities_each):
        """Compose the `misc` output node"""

        out_dict = {}
        gather_quantities(quantities_each, self.user_config.file_mapping['vasprun.xml'], out_dict, MISC_QUANTITIES)
        gather_quantities(quantities_each, self.user_config.file_mapping['OUTCAR'], out_dict, MISC_QUANTITIES)
        gather_quantities(quantities_each, self.user_config.file_mapping['vasp_output'], out_dict, MISC_QUANTITIES)
        return orm.Dict(dict=out_dict)

    def _compose_structure(self, quantities_each):
        """Compose the `structure` output node"""

        data = None
        if 'vasprun.xml' in quantities_each:
            data = quantities_each['vasprun.xml'].get('structure')
        if data is None:
            data = quantities_each.get('CONTCAR', {}).get('structure')
        if data is None:
            raise QuantityMissingError()

        node = orm.StructureData()
        node.set_cell(data['unitcell'])
        for site in data['sites']:
            node.append_atom(position=site['position'], symbols=site['symbol'], name=site['kind_name'])
        return node

    def _compose_arrays(self, quantities_each):
        """Generate the generic `arrays` output node"""
        array_quantities = ('projectors', 'born_charges', 'dielectrics', 'hessian', 'dynmat', 'energies')
        out_arrays = {}
        gather_quantities(quantities_each, 'vasprun.xml', out_arrays, array_quantities, flatten_dict=True)
        # Remove None values in the arrays
        out_arrays = {key: value for key, value in out_arrays.items() if value is not None}
        if out_arrays:
            return orm.ArrayData(out_arrays)
        return None

    def _compose_kpoints(self, quantities_each):
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
            return node
        raise QuantityMissingError('No valid kpoints data to use')

    def _compose_trajectory(self, quantities_each):
        """Compose the `trajectory` output"""

        node = orm.TrajectoryData()

        if 'vasprun.xml' in quantities_each:
            traj_data = quantities_each['vasprun.xml'].get('trajectory')
            if traj_data is None:
                return None
            for item in traj_data:
                for key, value in traj_data[item].items():
                    if key == 'symbols':
                        node.base.attributes.set(key, value)
                    else:
                        node.set_array(key, value)
            return node
        return None

    def _compose_bands(self, quantities_each):
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
            return node

    def _compose_dos(self, quantities_each):
        """Compose the `dos` node"""
        arrays_dict = {}
        if 'vasprun.xml' in quantities_each:
            gather_quantities(quantities_each, 'dos', arrays_dict, ['dos'], flatten_dict=True)
        if arrays_dict:
            node = orm.ArrayData(arrays_dict['dos'])
            return node

    def _check_vasp_errors(self, parser_notifications):  # pylint: disable=too-many-return-statements
        """
        Detect simple vasp execution problems and returns the exit_codes to be set
        """
        quantities = {}
        for key, value in self._quantities_each.items():
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
                    quantities,
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


def gather_quantities(quantities_each, namespace, dst, fields, flatten_dict=False):
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

    def __init__(self, notifications, parsed_quantities, inputs, exit_codes, critical_notifications):
        """
        Composed error codes based on the notifications

        Some of the errors need to have additional properties inspected before they can be emitted,
        as they might be trigged in a harmless way.

        To add new checkers, one needs to implement a property with the name of the error for this class and
        contains the code for checking. This property should return the exit_code or return None. The property
        is inspected if its name is in the list critical notifications.

        :param notification: The list of parsed notifications from the stream parser.
        :param parsed_quantities: The dictionary of parsed quantities.
        :param inputs: The dictionary of the input nodes.
        :param exit_codes: The dictionary of the exit codes from the parser.
        :param ignored: A list of critical notification that are allowed to have
        """
        self.notifications = notifications
        self.notifications_dict = {item['name']: item['message'] for item in self.notifications}
        self.parsed_quantities = parsed_quantities
        self.inputs = inputs
        self.exit_codes = exit_codes
        self.critical_notifications = critical_notifications

    def compose(self):
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
    def brmix(self):
        """Check if BRMIX should be emitted"""
        if 'brmix' not in self.notifications_dict:
            return None

        # If NELECT is set explicitly for the calculation then this is not an critical error
        if 'parameters' in self.inputs and 'nelect' in self.inputs['parameters'].get_dict():
            return None

        return self.exit_codes.ERROR_VASP_CRITICAL_ERROR.format(error_message=self.notifications_dict['brmix'])

    @property
    def edddav_zhegv(self):
        """Check if EDDDAV call to ZHEGV should be emitted. Sometimes it has converged."""
        if 'edddav_zhegv' not in self.notifications_dict:
            return None

        if self.parsed_quantities['run_status']['electronic_converged']:
            return None

        return self.exit_codes.ERROR_VASP_CRITICAL_ERROR.format(error_message=self.notifications_dict['edddav_zhegv'])

    @property
    def eddrmm_zhegv(self):
        """Check if EDDRMM call to ZHEGV should be emitted. Sometimes it has converged."""
        if 'eddrmm_zhegv' not in self.notifications_dict:
            return None

        if self.parsed_quantities['run_status']['electronic_converged']:
            return None

        return self.exit_codes.ERROR_VASP_CRITICAL_ERROR.format(error_message=self.notifications_dict['eddrmm_zhegv'])
