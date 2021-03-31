"""
VASP parser.

------------
The main driver routine for parsing VASP related files. The parser is very modular and
contains several modules:

- ``node_composer`` handles the quantity composition of nodes
- ``quantity`` the actual quantity to parse and what file parsers to use to obtain it
- ``settings`` general parser settings
- ``manager`` takes the quantity definitions and executes the actual parsing needed
"""
#encoding: utf-8
# pylint: disable=no-member
import traceback

from aiida.common.exceptions import NotExistent
from aiida_vasp.parsers.base import BaseParser
from aiida_vasp.parsers.quantity import ParsableQuantities
from aiida_vasp.parsers.settings import ParserSettings, ParserDefinitions
from aiida_vasp.parsers.node_composer import NodeComposer, get_node_composer_inputs

DEFAULT_OPTIONS = {
    'add_trajectory': False,
    'add_bands': False,
    'add_chgcar': False,
    'add_dos': False,
    'add_kpoints': False,
    'add_energies': False,
    'add_misc': True,
    'add_structure': False,
    'add_projectors': False,
    'add_born_charges': False,
    'add_dielectrics': False,
    'add_hessian': False,
    'add_dynmat': False,
    'add_wavecar': False,
    'add_forces': False,
    'add_stress': False,
    'add_site_magnetization': False,
}


class VaspParser(BaseParser):
    """
    Parses all VASP calculations.

    This particular class manages all the specific file parsers in
    aiida_vasp.parsers.file_parsers. The parser will check which quantities to parse
    and which nodes to add to the calculation based on the 'parser_settings' card in
    the 'settings' Dict of the corresponding VaspCalculation.

    Parser Settings usage:

    Parser settings can be passed through the input node `settings` as follows::

        settings = Dict(dict={
            'parser_settings': {
                ...
            }
        })

    Valid keys for `parser_settings` are:

    * `add_<quantity>`, where quantity is one of:

        'misc': (Default) Parameterdata node containing various quantities from OUTCAR and vasprun.xml.
        'structure':  (Default) StructureData node parsed from CONTCAR
        'bands':      Band structure node parsed from EIGENVAL.
        'dos':        ArrayData node containing the DOS parsed from DOSCAR.
        'kpoints':    KpointsData node parsed from IBZKPT.
        'wavecar':    FileData node containing the WAVECAR file.
        'chgcar':     FileData node containing the CHGCAR file.

    * `output_params`: A list of quantities, that should be added to the 'misc' node.

    * `file_parser_set`: String (DEFAULT = 'default').

        By this option the default set of FileParsers can be chosen. See settings.py
        for available options.

    Additional FileParsers can be added to the VaspParser by using

        VaspParser.add_file_parser(parser_name, parser_definition_dict),

    where the 'parser_definition_dict' should contain the 'parser_class' and the
    'is_critical' flag. Keep in mind adding an additional FileParser after 'parse_with_retrieved'
    is called, will only have an effect when parsing a second time.
    """

    def __init__(self, node):
        super(VaspParser, self).__init__(node)

        try:
            calc_settings = self.node.inputs.settings
        except NotExistent:
            calc_settings = None

        parser_settings = None
        if calc_settings:
            parser_settings = calc_settings.get_dict().get('parser_settings')

        self._definitions = ParserDefinitions()
        self._settings = ParserSettings(parser_settings, default_settings=DEFAULT_OPTIONS)
        self._parsable_quantities = ParsableQuantities(vasp_parser_logger=self.logger)
        self._file_parse_exit_codes = {}

    def add_parser_definition(self, filename, parser_dict):
        """Add the definition of a fileParser to self._definitions."""
        self._definitions.add_parser_definition(filename, parser_dict)

    def add_parsable_quantity(self, quantity_name, quantity_dict):
        """Add a single parsable quantity to the _parsable_quantities."""
        self._parsable_quantities.add_parsable_quantity(quantity_name, quantity_dict)

    def add_custom_node(self, node_name, node_dict):
        """Add a custom node to the settings."""
        self._settings.add_output_node(node_name, node_dict)

    def parse(self, **kwargs):
        """The function that triggers the parsing of a calculation."""

        self._file_parse_exit_codes = {}
        error_code = self._compose_retrieved_content(kwargs)
        if error_code is not None:
            return error_code

        for file_name, value_dict in self._definitions.parser_definitions.items():
            if file_name not in self._retrieved_content.keys() and value_dict['is_critical']:
                return self.exit_codes.ERROR_CRITICAL_MISSING_FILE

        self._parsable_quantities.setup(retrieved_filenames=self._retrieved_content.keys(),
                                        parser_definitions=self._definitions.parser_definitions,
                                        quantity_names_to_parse=self._settings.quantity_names_to_parse)

        # Parse the quantities from retrived files
        parsed_quantities, failed_to_parse_quantities = self._parse_quantities()

        # Store any exit codes returned in parser_warnings
        if self._file_parse_exit_codes:
            parsed_quantities['file_parser_warnings'] = self.parser_warnings

        # Compose the output nodes using the parsed quantities
        nodes_failed_to_create = self._compose_nodes(parsed_quantities)

        # Deal with missing quantities
        if failed_to_parse_quantities:
            return self.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY.format(quantity=', '.join(failed_to_parse_quantities))

        # Deal with missing node/nodes
        if nodes_failed_to_create:
            return self.exit_codes.ERROR_NOT_ABLE_TO_CREATE_NODE.format(nodes=', '.join(nodes_failed_to_create))

        # All quantities has been parsed, but there exit_codes reported from the parser
        # in this case, we return the code with the lowest status (hopefully the most severe)
        if self._file_parse_exit_codes:
            self._file_parse_exit_codes.sort(key=lambda x: x.status)
            return self._file_parse_exit_codes[0]

        return self.exit_codes.NO_ERROR

    def _parse_quantities(self):
        """
        This method dispatch the parsing to file parsers

        :returns: A tuple of parsed quantities dictionary and a list of quantities failed to obtain due to exceptions
        """
        parsed_quantities = {}
        # A dictionary for catching instantiated file parser objects
        file_parser_instances = {}
        failed_to_parse_quantities = []
        for quantity_key in self._parsable_quantities.quantity_keys_to_parse:
            file_name = self._parsable_quantities.quantity_keys_to_filenames[quantity_key]
            file_parser_cls = self._definitions.parser_definitions[file_name]['parser_class']

            # If a parse object has been instantiated, use it.
            if file_parser_cls in file_parser_instances:
                parser = file_parser_instances[file_parser_cls]
            else:
                try:
                    # The next line may except for ill-formated file
                    parser = file_parser_cls(settings=self._settings, exit_codes=self.exit_codes, file_path=self._get_file(file_name))
                except Exception:  # pylint: disable=broad-except
                    parser = None
                    failed_to_parse_quantities.append(quantity_key)
                    self.logger.warning('Cannot instantiate {}, exception {}:'.format(quantity_key, traceback.format_exc()))

                file_parser_instances[file_parser_cls] = parser

            # if the parser cannot be instantiated, add the quantity to a list of unavalaible ones
            if parser is None:
                failed_to_parse_quantities.append(quantity_key)
                continue

            # The next line may still except for ill-formated file - some parser load all data at
            # instantiation time, the others may not. See the `BaseFileParser.get_quantity`
            try:
                # The next line may still except for ill-formated file - some parser load all data at
                # instantiation time, the others may not
                parsed_quantity = parser.get_quantity(quantity_key)
            except Exception:  # pylint: disable=broad-except
                parsed_quantity = None
                failed_to_parse_quantities.append(quantity_key)
                self.logger.warning('Error parsing {} from {}, exception {}:'.format(quantity_key, parser, traceback.format_exc()))

            if parsed_quantity is not None:
                parsed_quantities[quantity_key] = parsed_quantity

            # Keep track of exit_code, if any
            if parser.exit_code and parser.exit_code.status != 0:
                self._file_parse_exit_codes[str(file_parser_cls)] = parser.exit_code

        return parsed_quantities, failed_to_parse_quantities

    def _compose_nodes(self, parsed_quantities):
        """
        Compose the nodes according to parsed quantities

        :returns: A list of link_names for the nodes that failed to compose
        """
        nodes_failed_to_create = []

        # Get the dictionary of equivalent quantities, and add a special quantity "parser_warnings"
        equivalent_quantity_keys = dict(self._parsable_quantities.equivalent_quantity_keys)

        for node_name, node_dict in self._settings.output_nodes_dict.items():
            inputs = get_node_composer_inputs(equivalent_quantity_keys, parsed_quantities, node_dict['quantities'])

            if node_name == 'misc':
                inputs['file_parser_warnings'] = parsed_quantities.get('file_parser_warnings')

            # If the input is empty, we skip creating the node as it is bound to fail
            if not inputs:
                nodes_failed_to_create.append(node_name)
                continue

            # Guard the parsing in case of errors
            try:
                aiida_node = NodeComposer.compose(node_dict['type'], inputs)
            except Exception:  # pylint: disable=broad-except
                nodes_failed_to_create.append(node_dict['link_name'])
                aiida_node = None
                self.logger.warning('Error creating output {} with type {}, exception: {}'.format(node_dict['link_name'], node_dict['type'],
                                                                                                  traceback.format_exc()))

            if aiida_node is not None:
                self.out(node_dict['link_name'], aiida_node)
        return nodes_failed_to_create

    @property
    def parser_warnings(self):
        """
        Compose a list of parser warnings as returned by individual file parsers
        """
        warnings = {}
        for key, exit_code in self._file_parse_exit_codes.items():
            warnings[key] = {
                'status': exit_code.status,
                'message': exit_code.message,
            }
        return warnings
