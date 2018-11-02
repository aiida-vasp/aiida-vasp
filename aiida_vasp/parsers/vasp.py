#encoding: utf-8

# pylint: disable=no-member
# Reason: pylint erroneously complains about non existing member 'get_quantity', which will be set in __init__.
"""AiiDA Parser for a aiida_vasp.VaspCalculation"""

from aiida_vasp.parsers.base import BaseParser
from aiida_vasp.parsers.quantity import ParsableQuantities, NODES
from aiida_vasp.parsers.parser_manager import ParserManager
from aiida_vasp.parsers.parser_settings import ParserSettings
from aiida_vasp.utils.delegates import Delegate

DEFAULT_OPTIONS = {
    'add_trajectory': False,
    'add_bands': False,
    'add_chgcar': False,
    'add_dos': False,
    'add_kpoints': False,
    'add_energies': False,
    'add_parameters': True,
    'add_structure': False,
    'add_projectors': False,
    'add_born_charges': False,
    'add_dielectrics': False,
    'add_hessian': False,
    'add_dynmat': False,
    'add_wavecar': False,
    'add_forces': False,
    'file_parser_set': 'default',
}


class VaspParser(BaseParser):
    """
    Parses all Vasp calculations.

    This particular class manages all the specific file parsers in
    aiida_vasp.io. The parser will check which quantities to parse and which nodes to add
    to the calculation based on the 'parser_settings' card in the 'settings' ParameterData of the
    corresponding VaspCalculation.

    Parser Settings usage:

    Parser settings can be passed through the input node `settings` as follows::

        settings = ParameterData(dict={
            'parser_settings': {
                ...
            }
        })

    Valid keys for `parser_settings` are:

    * `add_<quantity>`, where quantity is one of:

        'parameters': (Default) Parameterdata node containing various quantities from OUTCAR and vasprun.xml.
        'structure':  (Default) StructureData node parsed from CONTCAR
        'bands':      Band structure node parsed from EIGENVAL.
        'dos':        ArrayData node containing the DOS parsed from DOSCAR.
        'kpoints':    KpointsData node parsed from IBZKPT.
        'wavecar':    FileData node containing the WAVECAR file.
        'chgcar':     FileData node containing the CHGCAR file.

    * `output_params`: A list of quantities, that should be added to the 'output_parameters' node.

    * `file_parser_set`: String (DEFAULT = 'default').

        By this option the default set of FileParsers can be chosen. See parser_settings.py
        for available options.

    Additional FileParsers can be added to the VaspParser by using

        VaspParser.add_file_parser(parser_name, parser_definition_dict),

    where the 'parser_definition_dict' should contain the 'parser_class' and the
    'is_critical' flag. Keep in mind adding an additional FileParser after 'parse_with_retrieved'
    is called, will only have an effect when parsing a second time.
    """

    def __init__(self, calc):
        super(VaspParser, self).__init__(calc)

        # Initialise the 'get_quantity' delegate:
        setattr(self, 'get_quantity', Delegate())

        self.out_folder = None

        calc_settings = self._calc.get_inputs_dict().get('settings')
        settings = None
        if calc_settings:
            settings = calc_settings.get_dict().get('parser_settings')
        self.settings = ParserSettings(settings, DEFAULT_OPTIONS)

        self.quantities = ParsableQuantities(vasp_parser=self)
        self.parsers = ParserManager(vasp_parser=self)

        self._output_nodes = {}

        # this list is for bookkeeping, to check whether a quantity has been requested
        # twice during the parsing cycle.
        self._requested_quantities = []

    def add_file_parser(self, parser_name, parser_dict):
        """Add the definition of a fileParser to self.settings and self.parsers."""

        self.settings.parser_definitions[parser_name] = parser_dict
        self.parsers.add_file_parser(parser_name, parser_dict)

    def add_parsable_quantity(self, quantity_name, quantity_dict):
        """Add a single parsable quantity to the _parsable_quantities."""
        self.quantities.additional_quantities[quantity_name] = quantity_dict

    def parse_with_retrieved(self, retrieved):

        def missing_critical_file():
            for file_name, value_dict in self.settings.parser_definitions.items():
                if file_name not in self.out_folder.get_folder_list() and value_dict['is_critical']:
                    return True
            return False

        self.check_state()
        self.out_folder = self.get_folder(retrieved)

        if not self.out_folder:
            return self.result(success=False)

        if missing_critical_file():
            # A critical file i.e. OUTCAR does not exist. Abort parsing.
            return self.result(success=False)

        # Get the _quantities from the FileParsers.
        self.quantities.setup()

        # Set the quantities to parse list. Warnings will be issued if a quantity should be parsed and
        # the corresponding files do not exist.
        self.parsers.setup()
        quantities_to_parse = self.parsers.get_quantities_to_parse()

        # Parse all implemented quantities in the quantities_to_parse list.
        while quantities_to_parse:
            quantity = quantities_to_parse.pop(0)
            self._output_nodes.update(self.get_quantity(quantity))

        # Add output nodes if the corresponding data exists.
        for key, value in self._output_nodes.items():
            quantity = self.quantities.get_by_name(key)
            if not quantity.is_node:
                # This quantity does not represent a node, continue with the next one.
                continue

            if quantity.nodeName in self.settings.nodes and value is None:
                # One of the requested output nodes is None, parsing has not been successful.
                return self.result(success=False)

            if value:
                self._set_node(quantity.nodeName, value)

        # Reset the 'get_quantity' delegate
        self.get_quantity.clear()

        return self.result(success=True)

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
            # The quantity is not in the output_nodes. Try to parse it
            self._output_nodes.update(self.get_quantity(quantity))

        # parsing the quantity without requesting it a second time was successful, remove it from requested_quantities.
        self._requested_quantities.remove(quantity)

        # since the quantity has already been parsed now as an input, we don't have to parse it a second time later.
        self.parsers.remove(quantity)

        return self._output_nodes.get(quantity)

    def _set_node(self, node_name, node):
        """Wrapper for self.add_node, checking whether the Node is None and using the correct linkname"""

        if node is not None:
            self.add_node(NODES[node_name], node)
