#encoding: utf-8

# pylint: disable=no-member
# Reason: pylint erroneously complains about non existing member 'get_quantity', which will be set in __init__.
"""AiiDA Parser for a aiida_vasp.VaspCalculation"""

from aiida_vasp.parsers.base import BaseParser
from aiida_vasp.parsers.file_parser_definitions import get_file_parser_set
from aiida_vasp.utils.delegates import Delegate
from aiida_vasp.utils.extended_dicts import DictWithAttributes

LINKNAME_DICT = {
    'parameters': 'output_parameters',
    'kpoints': 'output_kpoints',
    'structure': 'output_structure',
    'trajectory': 'output_trajectory',
    'bands': 'output_bands',
    'dos': 'output_dos',
    'energies': 'output_energies',
    'projectors': 'output_projectors',
    'born_charges': 'output_born_charges',
    'dielectrics': 'output_dielectrics',
    'hessian': 'output_hessian',
    'dynmat': 'output_dynmat',
    'chgcar': 'chgcar',
    'wavecar': 'wavecar',
}

ALLOWED_TYPES = {
    'parameters': 'parameters',
    'kpoints': 'array.kpoints',
    'structure': 'structure',
    'trajectory': 'array.trajectory',
    'bands': 'array.bands',
    'dos': 'array',
    'energies': 'array',
    'projectors': 'array',
    'born_charges': 'array',
    'dielectrics': 'array',
    'hessian': 'array',
    'dynmat': 'array'
}

DEFAULT_OPTIONS = {
    'add_trajectory': False,
    'add_bands': False,
    'add_chgcar': False,
    'add_dos': False,
    'add_kpoints': True,
    'add_energies': True,
    'add_parameters': True,
    'add_structure': True,
    'add_projectors': False,
    'add_born_charges': False,
    'add_dielectrics': False,
    'add_hessian': False,
    'add_dynmat': False,
    'add_wavecar': False,
    'file_parser_set': 'default',
}


class ParsableQuantity(DictWithAttributes):
    """Container class for parsable quantities."""

    def __init__(self, name, init, files_list):
        # assign default values for optional attributes
        self.parsers = []
        self.alternatives = []

        super(ParsableQuantity, self).__init__(init)
        self.name = name

        # Check whether all files required for parsing this quantity have been retrieved and store it.
        missing_files = self.has_all(files_list)
        if missing_files is None:
            self.has_files = False
            missing_files = []
        else:
            self.has_files = not missing_files
        self.missing_files = missing_files

        # Check whether this quantity represents a node
        self.is_node = self.nodeName in LINKNAME_DICT

    def has_all(self, available_items):
        """Check whether all items are in item_list."""
        missing_items = []
        if not self.parsers:
            return None
        if not available_items:
            return self.parsers
        for item in self.parsers:
            if item not in available_items:
                missing_items.append(item)
        return missing_items


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

    * `file_parser_set`: String (DEFAULT = 'default').

        By this option the default set of FileParsers can be chosen. See file_parser_definitions.py
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

        self._parsable_quantities = {}

        self._settings = DEFAULT_OPTIONS
        calc_settings = self._calc.get_inputs_dict().get('settings')
        if calc_settings:
            self._settings.update(calc_settings.get_dict().get('parser_settings', DEFAULT_OPTIONS))

        file_parser_set = get_file_parser_set(self._settings['file_parser_set'])
        if file_parser_set is None:
            self.logger.warning('The FileParser set: {file_parser_set} has been '
                                'requested by setting `file_parser_set: {file_parser_set}`.'
                                'However it does not exist. Parsing will continue using '
                                'the `default` set of FileParser. Please check the '
                                '`parsers/file_parser_definitions.py` or the documentation '
                                'for available options.'.format(file_parser_set=self._settings['file_parser_set']))
            file_parser_set = get_file_parser_set()

        self._parsers = {}
        for key, value in file_parser_set.items():
            self.add_file_parser(key, value)

        self._quantities_to_parse = []
        self._output_nodes = {}

        # this list is for bookkeeping, to check whether a quantity has been requested
        # twice during the parsing cycle.
        self._requested_quantities = []

    def add_file_parser(self, parser_name, parser_dict):
        """
        Add the definition of a fileParser to self._parsers.

        The required tags for the parser_dict can be found in PARSABLE_FILES. The FileParser must inherit
        from BaseFileParser and it will replace another previously defined fileParser with the same name.
        """

        self._parsers[parser_name] = DictWithAttributes(parser_dict)

    def add_parsable_quantity(self, quantity_name, quantity_dict, retrieved_files=None):
        """Add a single parsable quantity to the _parsable_quantities."""

        self._parsable_quantities[quantity_name] = ParsableQuantity(quantity_name, quantity_dict, retrieved_files)

    def add_quantity_to_parse(self, quantity_to_add):
        """Check, whether a quantity or it's alternatives can be added."""
        for item in [quantity_to_add] + self._parsable_quantities[quantity_to_add].alternatives:
            if self._parsable_quantities[item].is_parsable:
                self._quantities_to_parse.append(item)
                return True
        return False

    def parse_with_retrieved(self, retrieved):

        def missing_critical_file():
            for file_name, value_dict in self._parsers.items():
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

        # Get the _parsable_quantities from the FileParsers.
        self._set_parsable_quantities()

        # Set the quantities to parse list. Warnings will be issued if a quantity should be parsed and
        # the corresponding files do not exist.
        self._check_and_validate_settings()

        self._set_file_parsers()

        # Parse all implemented quantities in the quantities_to_parse list.
        print self._quantities_to_parse
        while self._quantities_to_parse:
            quantity = self._quantities_to_parse.pop(0)
            balla = self.get_quantity(quantity, self._settings)
            if quantity == 'parameters':
                print self._settings
                print self._output_nodes
                print balla['parameters'].get_dict()
            self._output_nodes.update(self.get_quantity(quantity, self._settings))

        # Add output nodes if the corresponding data exists.
        for key, value in self._output_nodes.items():

            if not self._parsable_quantities[key].is_node:
                # This quantity does not represent a node, continue with the next one.
                continue

            if value:
                self._set_node(self._parsable_quantities[key].nodeName, value)

        return self.result(success=True)

    def _check_uniqueness_add_parsable(self):
        """Check uniqueness and add parsable quantities."""

        self._parsable_quantities = {}
        # Gather all parsable items as defined in the file parsers.
        for filename, value in self._parsers.items():
            # Initialise the instance of this FileParser to None.
            value.parser = None
            for quantity, quantity_dict in value['parser_class'].PARSABLE_ITEMS.items():
                if quantity in self._parsable_quantities:
                    # Check uniqueness
                    raise RuntimeError('The quantity {quantity} defined in {filename} has been '
                                       'defined by two FileParser classes. Quantity names must '
                                       'be unique. If both quantities are equivalent, define one '
                                       'as an alternative for the other.'.format(quantity=quantity, filename=filename))
                # Create quantity objects.
                self.add_parsable_quantity(quantity, quantity_dict, self.out_folder.get_folder_list())

    def _check_consitency_and_alternatives(self):
        """Check the consistency and alternatives."""

        import copy

        # Make a local copy of parsable_quantities, because during the next step
        # dummy quantities for missing quantities might be added.
        parsable_quantities = copy.deepcopy(self._parsable_quantities.keys())

        # Check for every quantity, whether all of the prerequisites are parsable, and whether
        # they are alternatives for another quantity.
        for quantity in parsable_quantities:
            value = self._parsable_quantities[quantity]
            is_parsable = True
            for prereq in value.prerequisites:
                if not self._parsable_quantities[prereq].has_files:
                    is_parsable = False
            value.is_parsable = is_parsable and value.has_files
            # Add this quantity to the list of alternatives of another quantity.
            if value.is_alternative is not None:
                if value.is_alternative not in self._parsable_quantities:
                    # The quantity which this quantity is an alternative to is not in _parsable_quantities.
                    # Add an unparsable dummy quantity for it.
                    self.add_parsable_quantity(value.is_alternative, {
                        'alternatives': [],
                        'nodeName': self._parsable_quantities[quantity].nodeName
                    })

                if quantity not in self._parsable_quantities[value.is_alternative].alternatives:
                    self._parsable_quantities[value.is_alternative].alternatives.append(quantity)

    def _set_parsable_quantities(self):
        """Set the parsable_quantities dictionary based on parsable_items obtained from the FileParsers."""

        # check uniqueness and add parsable quantities
        self._check_uniqueness_add_parsable()

        # check consistency, that the quantity is parsable and
        # alternatives
        self._check_consitency_and_alternatives()

    def _check_and_validate_settings(self):
        """Check the settings and set which files should be parsed based on the input."""

        self._quantities_to_parse = []
        for key, value in self._settings.items():
            if not key.startswith('add_'):
                # only keys starting with 'add_' will change the behaviour of the parser so get the next one.
                continue
            if not value:
                # The quantity should not be added, so the corresponding files do not have to be parsed.
                continue
            quantity = key[4:]
            if quantity not in self._parsable_quantities:
                self.logger.warning('{quantity} has been requested by setting '
                                    'add_{quantity}. However it has not been implemented. '
                                    'Please check the docstrings in aiida_vasp.parsers.vasp.py '
                                    'for valid input.'.format(quantity=quantity))
                continue

            # Found a node, which should be added, add it to the quantities to parse.
            # if all files required for this quantity have been retrieved. If there are
            # alternatives for this quantity also try those.
            success = self.add_quantity_to_parse(quantity)

            if not success:
                # Neither the quantity nor it's alternatives could be added to the quantities_to_parse.
                # Gather a list of all the missing files and issue a warning.
                missing_files = []
                for quant in [quantity] + self._parsable_quantities[quantity].alternatives:
                    for missing_file in self._parsable_quantities[quant].missing_files:
                        missing_files.append(missing_file)

                missing_files = ", ".join(missing_files)
                self.logger.warning('{quantity} has been requested, however the '
                                    'following files required for parsing have not been '
                                    'retrieved: {missing_files}.'.format(quantity=quantity, missing_files=missing_files))

    def _set_file_parsers(self):
        """
        Set the specific file parsers for OUTCAR, DOSCAR, EIGENVAL and vasprun.xml.

        Return False if a critical file is missing, which will abort the parsing.
        """

        for quantity in self._quantities_to_parse:
            for filename in self._parsable_quantities[quantity]['parsers']:
                if self._parsers[filename].parser is not None:
                    # This parser has already been checked, i.e. take the first
                    # available in the list that can be parsed (i.e. file exists)
                    continue
                file_to_parse = self.get_file(filename)
                self._parsers[filename].parser = self._parsers[filename]['parser_class'](self, file_path=file_to_parse)

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
            self._output_nodes.update(self.get_quantity(quantity, self._settings))

        # parsing the quantity without requesting it a second time was successful, remove it from requested_quantities.
        self._requested_quantities.remove(quantity)

        # since the quantity has already been parsed now as an input, we don't have to parse it a second time later.
        if quantity in self._quantities_to_parse:
            self._quantities_to_parse.remove(quantity)

        return self._output_nodes.get(quantity)

    def _set_node(self, node_name, node):
        """Wrapper for self.add_node, checking whether the Node is None and using the correct linkname"""

        if node is not None:
            self.add_node(LINKNAME_DICT[node_name], node)
