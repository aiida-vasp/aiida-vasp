from __future__ import annotations

from typing import Any

from aiida import orm
from aiida.engine import ExitCode

from aiida_vasp.parsers.content_parsers import *

from .vasp import MissingFileError, NotificationComposer, QuantityMissingError, VaspParser, get_structure_node

DEFAULT_EXCLUDED_QUANTITIES = (
    'energies',
    'chgcar',
    'wavecar',
    'projectors',
    'charge_density',
    'magnetization_density',
    'elastic_moduli',
    'symmetries',
)

DEFAULT_EXCLUDED_NODE = tuple(['bands', 'dos', 'kpoints', 'trajectory'])

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
    'neb_data',
    'outcar_positions',
    'outcar_cell',
)


class NebParser(VaspParser):
    """
    Parser for handling NEB calculations.
    """

    def parse(self, **kwargs: Any) -> ExitCode | None:
        user_config = self._init_user_settings()
        # Clear the critical objects list as they do not apply for NEB
        user_config.critical_objects = []

        exit_code = self._get_quantities_to_parse()
        if exit_code is not None:
            return exit_code

        # Check the number of images
        nimages = self.get_num_images()
        self.neb_indices = [f'{i:02d}' for i in range(1, nimages + 1)]

        # Parse the files
        def parse_and_add(
            name: str,
            parser_cls: Any,
            index: Any,
            required: bool = True,
            open_mode: str = 'r',
            content_parser_settings: dict | None = None,
        ) -> None:
            """
            Parse the target file and add the result to the quantities_each dictionary
            For NEB calculations, the results are in the individual folders named with the image index.
            The parsed quantities are rested inside a dictionary with the image index as the key.

            """
            resolved_name = user_config.file_mapping[name]
            # The stdout of the images are in the individual folders with the name stdout
            if name == 'vasp_output' and int(index) > 1:
                resolved_name = 'stdout'
            fpath = index + '/' + resolved_name
            # With the exception that the stdout of the first image is in the base calculation folder
            if name == 'vasp_output' and int(index) == 1:
                fpath = 'stdout'
            try:
                with self.retrieved.open(fpath, open_mode) as handler:
                    try:
                        parser: BaseFileParser = parser_cls(handler=handler, settings=content_parser_settings)
                    except Exception as error:
                        self.errored_parsers[name] = error
                        return
                    if parser.parser_notifications:
                        self.parser_notifications.update(parser.parser_notifications)
                    # Create the dictionary if it doesn't exist
                    if name not in self.quantities_each:
                        self.quantities_each[name] = {}
                    quantities, errored = parser.get_all_quantities()
                    # Save the parsed quantities in the dictionary with index as the key
                    for key, value in quantities.items():
                        if key not in self.quantities_each[name]:
                            self.quantities_each[name][key] = {}
                        self.quantities_each[name][key][index] = value
                    self.errored_quantities.update(errored)
            except FileNotFoundError as error:
                raise MissingFileError(f'{fpath} is missing in the retrieved folder.') from error

        for index in self.neb_indices:
            # We do not parse vasprun.xml as it is not always present in the sub folders
            # TODO - add XDATACAR support
            parse_and_add('OUTCAR', VtstNebOutcarParser, index, required=True)
            parse_and_add('vasp_output', StreamParser, index, required=True)
            parse_and_add('CONTCAR', PoscarParser, index, required=True)

        exit_code = self._post_process_quantities()
        if exit_code is not None:
            return exit_code

        return self._create_outputs()

    def get_num_images(self) -> int:
        """
        Return the number of images
        """
        try:
            nimages = self.node.inputs.parameters['images']
        except KeyError as no_images:
            raise ValueError('No `images` key defined in inputs - this is really an NEB calculation?') from no_images
        return nimages

    def _create_outputs(self) -> ExitCode | None:
        """Create the outputs"""
        # Create the outputs
        self._failed_to_compose = {}
        for name in ['misc', 'structure', 'trajectory', 'arrays']:
            if name in self.nodes_to_exclude:
                continue
            item = None
            try:
                item = getattr(self, '_compose_' + name)(self.quantities_each)
            except (QuantityMissingError, KeyError, ValueError, TypeError) as error:
                self._failed_to_compose[name] = error
                self.logger.warning(f'Failed to compose {name} node: {error}')
                continue
            # If a Node is returned, add it to the outputs
            if isinstance(item, orm.Data):
                self.out(name, item)
            # If a dictionary is returned, it is a namespace output
            if isinstance(item, dict):
                for key, value in item.items():
                    self.out(name + '.' + key, value)
        if (
            any(name in self.user_config.include_node for name in self._failed_to_compose)
            and self.user_config.check_completeness is True
        ):
            return self.exit_codes.ERROR_NOT_ABLE_TO_CREATE_NODE.format(nodes=', '.join(self._failed_to_compose.keys()))
        # Check for errors
        if self.user_config.check_errors is True:
            error = self._check_vasp_errors(self.parser_notifications)
            return error

    def _compose_structure(self, quantities_each: dict[str, Any]) -> dict[str, orm.StructureData]:
        """Compose the `structure` output nodes"""

        data = quantities_each['CONTCAR'].get('structure')
        if len(data) == 0:
            raise QuantityMissingError()

        output = {}
        for index in self.neb_indices:
            output['image_' + index] = get_structure_node(data[index])
        return output

    # At the moment there is no array output for NEB calculations
    # def _compose_arrays(self, quantities_each):
    #     """Generate the generic `arrays` output node"""
    #     array_quantities = ('energies',)
    #     out_arrays = {}
    #     gather_quantities_neb(quantities_each, 'OUTCAR', out_arrays, array_quantities)
    #     # Remove None values in the arrays
    #     out_arrays = {key: value for key, value in out_arrays.items() if value is not None}
    #     if out_arrays:
    #         return orm.ArrayData(out_arrays)
    #     return None

    def _compose_trajectory(self, quantities_each: dict[str, Any]) -> dict[str, orm.TrajectoryData] | None:
        """Compose the `trajectory` output node"""
        output = {}
        for index in self.neb_indices:
            node = orm.TrajectoryData()
            if 'vasprun.xml' in quantities_each:
                traj_data = quantities_each['vasprun.xml'].get('trajectory')
                if traj_data is None:
                    return None
                for key, value in traj_data.items():
                    if key == 'symbols':
                        node.base.attributes.set(key, value)
                    else:
                        node.set_array(key, value)
                output['image_' + index] = node
        return output

    def _compose_misc(self, quantities_each: dict[str, Any]) -> orm.Dict:
        """Compose the `misc` output node"""

        out_dict = {}
        gather_quantities_neb(quantities_each, self.user_config.file_mapping['vasprun.xml'], out_dict, MISC_QUANTITIES)
        gather_quantities_neb(quantities_each, self.user_config.file_mapping['OUTCAR'], out_dict, MISC_QUANTITIES)
        gather_quantities_neb(quantities_each, self.user_config.file_mapping['vasp_output'], out_dict, MISC_QUANTITIES)
        return orm.Dict(dict=out_dict)

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
        for index in self.neb_indices:
            run_status = quantities['run_status'][index]

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
                self.logger.warning(
                    'The ionic relaxation is not converged, but the calculation is treated as successful.'
                )

        # Check for the existence of critical warnings
        if 'notifications' in quantities:
            notifications = quantities['notifications']['01']
            ignore_all = self.user_config.ignore_notification_errors
            if not ignore_all:
                composer = NotificationComposer(
                    notifications,
                    quantities['run_status']['01'],
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


def gather_quantities_neb(
    quantities_each: dict[str, Any], namespace: str, dst: dict[str, Any], fields: list[str]
) -> None:
    """
    Gather quantities and put them into the target dictionary
    """
    for key, value in quantities_each.get(namespace, {}).items():
        if key in fields:
            dst[key] = value
