"""
Parser for NEB calculations using VASP compiled with VTST
"""
import os
import traceback
from pathlib import Path

from aiida.common.exceptions import NotExistent
from aiida_vasp.parsers.settings import ParserSettings, ParserDefinitions
from aiida_vasp.parsers.node_composer import NodeComposer
from aiida_vasp.parsers.quantity import ParsableQuantities
from aiida_vasp.parsers.vasp import VaspParser, NotificationComposer

# pylint: disable=no-member
# pylint: disable=logging-fstring-interpolation

NEB_NODES = {
    'neb_misc': {
        'link_name': 'neb_misc',
        'type': 'dict',
        'quantities': ['neb_data']
    },
    'misc': {
        'link_name': 'misc',
        'type': 'dict',
        'quantities': [
            'notifications',
            'run_stats',
            'run_status',
        ]
    },
    'kpoints': {
        'link_name': 'kpoints',
        'type': 'array.kpoints',
        'quantities': ['kpoints-kpoints'],
    },
    'structure': {
        'link_name': 'structure',
        'type': 'structure',
        'quantities': ['structure'],
    },
    'chgcar': {
        'link_name': 'chgcar',
        'type': 'vasp.chargedensity',
        'quantities': ['chgcar'],
    },
    'wavecar': {
        'link_name': 'wavecar',
        'type': 'vasp.wavefun',
        'quantities': ['wavecar'],
    },
    'site_magnetization': {
        'link_name': 'site_magnetization',
        'type': 'dict',
        'quantities': ['site_magnetization'],
    },
    'image_forces': {
        'link_name': 'image_forces',
        'type': 'array',
        'quantities': ['forces']
    }
}

DEFAULT_SETTINGS = {
    'add_bands': False,
    'add_chgcar': False,
    'add_dos': False,
    'add_kpoints': False,
    'add_misc': True,
    'add_neb_misc': True,
    'add_structure': True,
    'add_wavecar': False,
    'add_site_magnetization': False,
    'add_image_forces': False,
    'critical_notifications': {
        'add_brmix': True,
        'add_cnormn': True,
        'add_denmp': True,
        'add_dentet': True,
        'add_edddav_zhegv': True,
        'add_eddrmm_zhegv': True,
        'add_edwav': True,
        'add_fexcp': True,
        'add_fock_acc': True,
        'add_non_collinear': True,
        'add_not_hermitian': True,
        #add_psmaxn': True,
        'add_pzstein': True,
        'add_real_optlay': True,
        'add_rhosyg': True,
        'add_rspher': True,
        'add_set_indpw_full': True,
        'add_sgrcon': True,
        'add_no_potimm': True,
        'add_magmom': True,
    }
}

_VASP_OUTPUT = 'stdout'


class NEBNodeComposer(NodeComposer):
    """
    NodeComposer for NEB

    Some quantities are composed at a per-image basis
    """
    COMBINED_NODES = ['neb_misc', 'image_forces']
    COMBINED_QUANTITY = ['neb_data', 'outcar-forces']

    @staticmethod
    def get_image_str(image_idx):
        return f'{image_idx:02d}'

    def get_num_images(self):
        """Infer the number of image from the keys of the quantities dictionary"""
        indices = []
        for key in self._quantities.keys():
            try:
                value = int(key)
            except ValueError:
                continue
            indices.append(value)
        return max(indices)

    def compose_nodes(self):  # pylint: disable=too-many-branches
        """
        Compose the nodes as required.

        The major difference compared to the standard calculations is that NEB
        calculations have different images and most data are parsed at a per image
        basis. For example, each image would have its own output structure, bands, and
        kpoints.

        However, some output includes the data from images, and needs to be handed separately.
        """

        # Compose nodes for each image and buildthe combined quantity dictionary
        combined_quantities = {}
        for image_idx in range(1, self.get_num_images() + 1):
            quantity_dict = {}
            image_name = self.get_image_str(image_idx)
            for key, value in self._quantities[image_name].items():
                if key in self.COMBINED_QUANTITY:
                    # Combined excluded data into a single dictionary
                    # This gives something like {'neb_data': {'01': {...dict for image 1...}},
                    #                            'image_forces': {'image_forces_01': {...dict for image 1...}}}
                    # The prefix for forces is needed so that the composer will name each array like 'image_forces_xx'
                    if key == 'outcar-forces':
                        sub_key = 'forces_' + image_name
                    else:
                        sub_key = image_name
                    if key in combined_quantities:
                        combined_quantities[key][sub_key] = value
                    else:
                        combined_quantities[key] = {sub_key: value}
                else:
                    quantity_dict[key] = value

            self._compose_nodes_for_image(image_idx)

        # Update with combined quantities in the top level
        self._quantities.update(combined_quantities)

        # Deal with the combined data
        for node_name, node_dict in self._nodes.items():

            # Deal with only the nodes containing combined data
            if node_name not in self.COMBINED_NODES:
                continue

            inputs = self._set_input_quantities(node_dict['quantities'])
            # If the input is empty, we skip creating the node as it is bound to fail
            if not inputs:
                self._failed_to_create.append(node_name)
                self._logger.warning(f'Creating node {node_dict["link_name"]} of type {node_dict["type"]} failed. '
                                     'No parsed data available.')
                continue
            exception = None
            # Guard the parsing in case of errors
            try:
                node = self.compose_node(node_dict['type'], inputs)
            except Exception:  # pylint: disable=broad-except
                node = None
                exception = traceback.format_exc()

            if node is not None:
                self._created[node_dict['link_name']] = node
            else:
                self._logger.warning(f'Creating node {node_dict["link_name"]} of type {node_dict["type"]} failed, '
                                     f'exception: {exception}')
                self._failed_to_create.append(node_dict['link_name'])

    def _compose_nodes_for_image(self, image_idx):
        """
        Compose the nodes according to parsed quantities

        :returns: A list of link_names for the nodes that failed to compose
        """
        nodes_failed_to_create = []

        for node_name, node_dict in self._nodes.items():
            if node_name in self.COMBINED_NODES:
                continue

            inputs = self._set_input_quantities(node_dict['quantities'], image_idx)
            # If the input is empty, we skip creating the node as it is bound to fail
            if not inputs:
                self._failed_to_create.append(node_name)
                self._logger.warning(
                    f'Creating node {node_dict["link_name"]} of type {node_dict["type"]} failed for image {image_idx:02d}. '
                    'No parsed data available.')
                continue
            exception = None
            # Guard the parsing in case of errors
            try:
                node = self.compose_node(node_dict['type'], inputs)
            except Exception:  # pylint: disable=broad-except
                node = None
                exception = traceback.format_exc()

            link_name = node_dict['link_name'] + f'.image_{image_idx:02d}'
            if node is not None:
                # Suffix the output name with image id
                self._created[link_name] = node
            else:
                self._logger.warning(f'Creating node {link_name} of type {node_dict["type"]} failed, ' f'exception: {exception}')
                self._failed_to_create.append(link_name)

        return nodes_failed_to_create

    def _set_input_quantities(self, node_quantities, image_idx=None):  # pylint: disable=arguments-differ
        """Set the necessary input quantities for the node."""
        inputs = {}
        # Iterate over the quantities that is requested for this node
        for quantity in node_quantities:
            # Find this quantity's equivalent quantities
            equivalent_quantities = self._find_equivalent_quantities(quantity)
            if equivalent_quantities is not None:
                # Check if these are parsed and pick the first one if multiple exists
                # Get a dictionary of quantities to be searched from
                quantity_dict = self._quantities if image_idx is None else self._quantities[self.get_image_str(image_idx)]
                for equivalent_quantity in equivalent_quantities:
                    if equivalent_quantity in quantity_dict:
                        # Make sure we strip prefixes as the quantities can contain
                        # multiple equivalent keys, relevant only up to now.
                        new_key = equivalent_quantity
                        if '-' in equivalent_quantity:
                            new_key = equivalent_quantity.split('-')[1]
                        inputs[new_key] = quantity_dict[equivalent_quantity]
                        break
        return inputs


class NEBSettings(ParserSettings):
    """
    Settings for NEB calculations
    """

    NODES = NEB_NODES


class VtstNebParser(VaspParser):
    """
    Parser for parsing NEB calculations performed with VASP compiled with VTST tools.

    The major difference compared with standard VASP calculations is that the output files are placed
    in subfolders. With the only exception being `vasprun.xml` it is not clear what image this file is
    for.
    """
    COMPOSER_CLASS = NEBNodeComposer

    def __init__(self, node):
        super(VtstNebParser, self).__init__(node)
        try:
            calc_settings = self.node.inputs.settings
        except NotExistent:
            calc_settings = None

        parser_settings = None
        if calc_settings:
            parser_settings = calc_settings.get_dict().get('parser_settings')

        self._settings = NEBSettings(parser_settings, default_settings=DEFAULT_SETTINGS, vasp_parser_logger=self.logger)
        self._definitions = ParserDefinitions(object_parser_set='neb')
        self._parsable_quantities = NEBParsableQuantities(vasp_parser_logger=self.logger)

    def get_num_images(self):
        """
        Return the number of images
        """
        try:
            nimages = self.node.inputs.parameters['images']
        except KeyError:
            raise ValueError('No `images` key defined in inputs - this is really an NEB calculation?')
        return nimages

    def _setup_parsable(self):
        """Setup the parable quantities. For NEB calculations we collpase the folder structure"""
        filenames = {Path(fname).name for fname in self._retrieved_content}
        self._parsable_quantities.setup(retrieved_content=list(filenames),
                                        parser_definitions=self._definitions.parser_definitions,
                                        quantity_names_to_parse=self._settings.quantity_names_to_parse)

    def _parse_quantities(self):
        """
        Parse the quantities. This has to be done for each image

        Returns:
            a dictionary with keys like: '01', '02'... and values being the parsed quantities for each image
        """
        nimages = self.get_num_images()

        per_image_quantities = {}
        #per_image_failed_quantities = {}
        failed_quantities = []

        for image_idx in range(1, nimages + 1):
            quantities, failed = self._parse_quantities_for_image(image_idx)
            per_image_quantities[f'{image_idx:02d}'] = quantities
            #per_image_failed_quantities[f'{image_idx:02d}'] = failed
            failed_quantities.extend([f'image_{image_idx:02d}_{name}' for name in failed])

        return per_image_quantities, failed_quantities

    # Override super class methods
    def _parse_quantities_for_image(self, image_idx):
        """
        This method dispatch the parsing to file parsers

        :returns: A tuple of parsed quantities dictionary and a list of quantities failed to obtain due to exceptions
        """
        parsed_quantities = {}
        # A dictionary for catching instantiated file parser objects
        file_parser_instances = {}
        failed_to_parse_quantities = []
        for quantity_key in self._parsable_quantities.quantity_keys_to_parse:
            name = self._parsable_quantities.quantity_keys_to_content[quantity_key]

            # Skip vasprun.xml file that does not exists for each image as of vasp 5.4.4
            if name == 'vasprun.xml':
                continue

            # Full path of the file, including the image folder
            content_path = f'{image_idx:02d}/' + name

            # Special case - for the stdout of the first image, we have to parse from the root stdout
            if image_idx == 1 and name == _VASP_OUTPUT:
                content_path = _VASP_OUTPUT

            object_parser_cls = self._definitions.parser_definitions[name]['parser_class']

            # If a parse object has been instantiated, use it.
            if object_parser_cls in file_parser_instances:
                parser = file_parser_instances[object_parser_cls]
            else:
                try:
                    # The next line may except for ill-formated file
                    with self._get_handler(content_path) as handler:
                        parser = object_parser_cls(settings=self._settings, handler=handler)
                except Exception:  # pylint: disable=broad-except
                    parser = None
                    failed_to_parse_quantities.append(quantity_key)
                    self.logger.warning('Cannot instantiate {} for image {}, exception {}:'.format(
                        object_parser_cls, image_idx, traceback.format_exc()))

                file_parser_instances[object_parser_cls] = parser

            # if the parser cannot be instantiated, add the quantity to a list of unavalaible ones
            if parser is None:
                failed_to_parse_quantities.append(quantity_key)
                continue

            # The next line may still except for ill-formated file - some parser load all data at
            # instantiation time, the others may not. See the `BaseFileParser.get_quantity`
            exception = None
            try:
                # The next line may still except for ill-formated file - some parser load all data at
                # instantiation time, the others may not
                parsed_quantity = parser.get_quantity(quantity_key)
            except Exception:  # pylint: disable=broad-except
                parsed_quantity = None
                exception = traceback.format_exc()

            if parsed_quantity is not None:
                parsed_quantities[quantity_key] = parsed_quantity
            else:
                self.logger.warning('Parsing {} from {} failed for image {}, exception: {}'.format(
                    quantity_key, parser, image_idx, exception))
                failed_to_parse_quantities.append(quantity_key)

        return parsed_quantities, failed_to_parse_quantities

    def _check_vasp_errors(self, quantities):
        """
        Detect simple vasp execution problems and returns the exit_codes to be set
        """

        # Check if some diagnosis information is missing
        neb_data_list = [image.get('neb_data') for image in quantities.values()]
        run_status_list = [image.get('run_status') for image in quantities.values()]

        if any(data is None for data in neb_data_list) or any(data is None for data in run_status_list):
            return self.exit_codes.ERROR_DIAGNOSIS_OUTPUTS_MISSING

        # Return errors related to execution and convergence problems.
        # Note that the order is important here - if a calculation is not finished, we cannot
        # comment on wether properties are converged are not.
        # Here we only check for the first frame
        if any(run_status['finished'] is False for run_status in run_status_list):
            return self.exit_codes.ERROR_DID_NOT_FINISH

        if run_status_list[0]['electronic_converged'] is False:
            return self.exit_codes.ERROR_ELECTRONIC_NOT_CONVERGED

        # Check the ionic convergence issues - the system is converged only if all image are "neb converged"
        if not all(per_image.get('neb_converged', False) for per_image in neb_data_list):
            if self._check_ionic_convergence:
                return self.exit_codes.ERROR_IONIC_NOT_CONVERGED
            self.logger.warning('The NEB calculation is not converged, but the calculation is treated as successful.')

        # Check for the existence of critical warnings. This is performed for all images.
        all_notifications = []
        for image in quantities.values():
            all_notifications.extend(image.get('notifications', []))

        ignore_all = self.parser_settings.get('ignore_all_errors', False)
        if not ignore_all:
            composer = NotificationComposer(all_notifications,
                                            quantities,
                                            self.node.inputs,
                                            self.exit_codes,
                                            parser_settings=self._settings)
            exit_code = composer.compose()
            if exit_code is not None:
                return exit_code

        return None


class NEBParsableQuantities(ParsableQuantities):
    """Quantity settings for NEB"""

    def _identify_missing_content(self, retrieved_content, names_in_parser_definitions):
        """Identify missing objects for quantities. Not there is no per-image checks"""
        _missing_content = {}
        # Check also sub folders
        filenames = [os.path.split(name)[-1] for name in retrieved_content]
        for quantity_key in self._quantity_items:
            name = self._quantity_keys_to_content[quantity_key]
            if name not in filenames or name not in names_in_parser_definitions:
                _missing_content[quantity_key] = name
        return _missing_content
