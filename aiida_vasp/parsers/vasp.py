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
from aiida_vasp.parsers.manager import ParserManager
from aiida_vasp.parsers.settings import ParserSettings
from aiida_vasp.parsers.node_composer import NodeComposer
from aiida_vasp.utils.delegates import Delegate

# defaults

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
    'file_parser_set': 'default',
    'store_energies_sc': False,
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

        By this option the default set of FileParsers can be chosen. See parser_settings.py
        for available options.

    Additional FileParsers can be added to the VaspParser by using

        VaspParser.add_file_parser(parser_name, parser_definition_dict),

    where the 'parser_definition_dict' should contain the 'parser_class' and the
    'is_critical' flag. Keep in mind adding an additional FileParser after 'parse_with_retrieved'
    is called, will only have an effect when parsing a second time.
    """

    def __init__(self, node):
        super(VaspParser, self).__init__(node)

        # Initialise the 'get_quantity' delegate:
        setattr(self, 'get_quantity', Delegate())

        try:
            calc_settings = self.node.inputs.settings
        except NotExistent:
            calc_settings = None

        parser_settings = None
        if calc_settings:
            parser_settings = calc_settings.get_dict().get('parser_settings')

        self.settings = ParserSettings(parser_settings, DEFAULT_OPTIONS)
        self.quantities = ParsableQuantities()
        self.parser_manager = ParserManager(node=self.node, vasp_parser_logger=self.logger)
        self._output_nodes = {}

        # this list is for bookkeeping, to check whether a quantity has been requested
        # twice during the parsing cycle.
        self._requested_quantities = []

    def add_file_parser(self, parser_name, parser_dict):
        """Add the definition of a fileParser to self.settings and self.parser_manager."""

        self.settings.parser_definitions[parser_name] = parser_dict
        self.parser_manager.add_file_parser(parser_name, parser_dict)

    def add_parsable_quantity(self, quantity_name, quantity_dict):
        """Add a single parsable quantity to the _parsable_quantities."""
        self.quantities.additional_quantities[quantity_name] = quantity_dict

    def add_custom_node(self, node_name, node_dict):
        """Add a custom node to the settings."""
        self.settings.add_output_node(node_name, node_dict)

    def parse(self, **kwargs):
        """The function that triggers the parsing of a calculation."""

        def missing_critical_file():
            for file_name, value_dict in self.settings.parser_definitions.items():
                if file_name not in self._retrieved_content.keys() and value_dict['is_critical']:
                    return True
            return False

        error_code = self.compose_retrieved_content(kwargs)

        if error_code is not None:
            return error_code
        if missing_critical_file():
            # A critical file is missing. Abort parsing
            # in case we do not find this or other files marked with is_critical
            return self.exit_codes.ERROR_CRITICAL_MISSING_FILE

        # Get the _quantities from the FileParsers.
        self.quantities.setup(retrieved_filenames=self._retrieved_content.keys(), parser_definitions=self.settings.parser_definitions)

        # Set the quantities to parse list. Warnings will be issued if a quantity should be parsed and
        # the corresponding files do not exist.
        self.parser_manager.setup(parser_definitions=self.settings.parser_definitions,
                                  quantities_to_parse=self.settings.quantities_to_parse,
                                  quantities=self.quantities)

        # Parse all implemented quantities in the quantities_to_parse list.
        quantity_name_to_file_name = {}
        for quantity_name in self.parser_manager.quantities_to_parse:
            file_name = self.quantities.get_by_name(quantity_name).file_name
            if quantity_name not in quantity_name_to_file_name:
                quantity_name_to_file_name[quantity_name] = file_name

        for quantity_name, file_name in quantity_name_to_file_name.items():
            file_to_parse = self.get_file(file_name)
            parser = self.parser_manager.parsers[file_name]['parser_class'](settings=self.settings,
                                                                            exit_codes=self.exit_codes,
                                                                            file_path=file_to_parse)
            parsed_data = parser.get_quantity(quantity_name)
            if parsed_data and parsed_data[quantity_name] is not None:
                self._output_nodes.update(parsed_data)
            self.exit_code = parser.exit_code

        # Assemble the nodes associated with the quantities
        node_assembler = NodeComposer(vasp_parser=self)
        for node_name, node_dict in self.settings.output_nodes_dict.items():
            node = node_assembler.compose(node_dict.type, node_dict.quantities)
            success = self._set_node(node_name, node)
            if not success:
                return self.exit_codes.ERROR_PARSING_FILE_FAILED

        # Reset the 'get_quantity' delegate
        self.get_quantity.clear()

        try:
            return self.exit_code
        except AttributeError:
            pass

        return self.exit_codes.NO_ERROR

    def get_inputs(self, quantity):
        """
        Return a quantity required as input for another quantity.

        This method will be called by the FileParsers in order to get a required input quantity
        from self._output_nodes. If the quantity is not in the dictionary the VaspParser will
        try to parse it. If a quantiy has been requested this way two times, parsing will be
        aborted because there is a cyclic dependency of the parsable items.
        """
        if quantity in self._requested_quantities:
            raise RuntimeError('{quantity} has been requested for parsing a second time. '
                               'There is probably a cycle in the prerequisites of the '
                               'parsable_items in the single FileParsers.'.format(quantity=quantity))

        # This is the first time this quantity has been requested, keep track of it.
        self._requested_quantities.append(quantity)
        if quantity not in self._output_nodes:
            # Did we parse an alternative
            for item in self.quantities.get_equivalent_quantities(quantity):
                if item.original_name in self._output_nodes:
                    return {quantity: self._output_nodes.get(item.original_name)}
            # The quantity is not in the output_nodes. Try to parse it
            self._output_nodes.update(self.get_quantity(quantity))

        # parsing the quantity without requesting it a second time was successful, remove it from requested_quantities.
        self._requested_quantities.remove(quantity)

        # since the quantity has already been parsed now as an input, we don't have to parse it a second time later.
        self.parser_manager.remove(quantity)

        return {quantity: self._output_nodes.get(quantity)}

    def _set_node(self, node_name, node):
        """Wrapper for self.add_node, checking whether the Node is None and using the correct linkname."""

        if node is None:
            return False
        self.out(self.settings.output_nodes_dict[node_name].link_name, node)
        return True
