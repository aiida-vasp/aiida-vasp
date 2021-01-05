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

        exit_code = None
        error_code = self._compose_retrieved_content(kwargs)
        if error_code is not None:
            return error_code

        for file_name, value_dict in self._definitions.parser_definitions.items():
            if file_name not in self._retrieved_content.keys() and value_dict['is_critical']:
                return self.exit_codes.ERROR_CRITICAL_MISSING_FILE

        self._parsable_quantities.setup(retrieved_filenames=self._retrieved_content.keys(),
                                        parser_definitions=self._definitions.parser_definitions,
                                        quantity_names_to_parse=self._settings.quantity_names_to_parse)

        parsed_quantities = {}
        for quantity_key in self._parsable_quantities.quantity_keys_to_parse:
            file_name = self._parsable_quantities.quantity_keys_to_filenames[quantity_key]
            file_parser_cls = self._definitions.parser_definitions[file_name]['parser_class']
            parser = file_parser_cls(settings=self._settings, exit_codes=self.exit_codes, file_path=self._get_file(file_name))
            parsed_quantity = parser.get_quantity(quantity_key)
            if parsed_quantity is not None:
                parsed_quantities[quantity_key] = parsed_quantity
            exit_code = parser.exit_code

        for _, node_dict in self._settings.output_nodes_dict.items():
            equivalent_quantity_keys = self._parsable_quantities.equivalent_quantity_keys
            inputs = get_node_composer_inputs(equivalent_quantity_keys, parsed_quantities, node_dict['quantities'])
            aiida_node = NodeComposer.compose(node_dict['type'], inputs)
            if aiida_node is None:
                return self.exit_codes.ERROR_PARSING_FILE_FAILED
            self.out(node_dict['link_name'], aiida_node)

        if exit_code is not None:
            return exit_code

        return self.exit_codes.NO_ERROR
