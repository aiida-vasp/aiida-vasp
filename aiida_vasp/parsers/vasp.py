#encoding: utf-8
"""AiiDA Parser for a aiida_vasp.VaspCalculation"""

from aiida_vasp.parsers.base import BaseParser
from aiida_vasp.parsers.file_parser_definitions import get_file_parser_set
from aiida_vasp.utils.delegates import delegate
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.array.bands import BandsData
from aiida.orm.data.array.trajectory import TrajectoryData
from aiida.orm.data.array import ArrayData
from aiida_vasp.utils.extended_dicts import DictWithAttributes

LINKNAME_DICT = {
    'parameters': 'output_parameters',
    'kpoints': 'output_kpoints',
    'structure': 'output_structure',
    'trajectory': 'output_trajectory',
    'trajectory_full': 'output_trajectory_arr',
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
    'parameters': ParameterData,
    'kpoints': KpointsData,
    'structure': StructureData,
    'trajectory': TrajectoryData,
    'trajectory_full': ArrayData,
    'bands': BandsData,
    'dos': ArrayData,
    'energies': ArrayData,
    'projectors': ArrayData,
    'born_charges': ArrayData,
    'dielectrics': ArrayData,
    'hessian': ArrayData,
    'dynmat': ArrayData
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

    QUANTITY_ATTR_DEFAULTS = {
        'alternatives': [],
        'parsers': [],
    }

    def __init__(self, name, init, files_list):
        super(ParsableQuantity, self).__init__(init, ParsableQuantity.QUANTITY_ATTR_DEFAULTS)
        self.name = name

        # Check whether all files required for parsing this quantity have been retrieved and store it.
        missing_files = self.has_all(files_list)
        if missing_files is None:
            self.has_files = False
            missing_files = []
        else:
            self.has_files = not missing_files
        self.missing_files = missing_files

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

        self.out_folder = None

        self._parsable_quantities = {}

        self._settings = DEFAULT_OPTIONS
        calc_settings = self._calc.get_inputs_dict().get('settings')
        if calc_settings:
            self._settings.update(calc_settings.get_dict().get('parser_settings', DEFAULT_OPTIONS))

        file_parser_set = get_file_parser_set(self._settings['file_parser_set'])
        if file_parser_set is None:
            self.logger.warning(' The {0} FileParser-set has been requested by setting `file_parser_set: {0}`,'.format(
                self._settings['file_parser_set']) + ' however it does not exist. Parsing will continue using the `default` set of' +
                                ' Please check the `parsers/file_parser_definitions.py`' + ' or the documentation for available options.')
            file_parser_set = get_file_parser_set()

        self._parsers = {}
        for key, value in file_parser_set.iteritems():
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

    def parse_with_retrieved(self, retrieved):

        def missing_critical_file():
            for file_name, value_dict in self._parsers.iteritems():
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
        print("quantities_to_parse:", self._quantities_to_parse)
        while self._quantities_to_parse:
            quantity = self._quantities_to_parse.pop(0)
            self._output_nodes.update(self.get_quantity(quantity, self._settings))
            
        # Add output nodes if the corresponding data exists.
        for key, value in self._output_nodes.iteritems():

            if not self._parsable_quantities[key].is_node:
                # This quantity does not represent a node, continue with the next one.
                continue

            if value:
                self._set_node(key, value)

        return self.result(success=True)

    def _set_parsable_quantities(self):
        """Set the parsable_quantities dictionary based on parsable_items obtained from the FileParsers."""

        import copy

        # Gather all parsable items as defined in the file parsers.
        for filename, value in self._parsers.iteritems():
            # initialise the instance of this FileParser to None.
            if filename in self._parsable_quantities:
                raise RuntimeError('The quantity {0} has been defined by two FileParser classes.'.format(filename) +
                                   ' Quantity names must be unique. If both quantities are equivalent, define one as' +
                                   ' an alternative for the other.')
            value.parser = None
            for quantity, quantity_dict in value['parser_class'].PARSABLE_ITEMS.iteritems():
                # Create quantity objects.
                self.add_parsable_quantity(quantity, quantity_dict, self.out_folder.get_folder_list())

        # make a local copy of parsable_quantities, because during the next step
        # dummy quantities for missing quantities might be added.
        parsable_quantities = copy.deepcopy(self._parsable_quantities.keys())

        # check for every quantity, whether all of the prerequisites are parsable, whether they
        # represent a node, and whether they are alternatives for another quantity.
        for quantity in parsable_quantities:
            value = self._parsable_quantities[quantity]
            is_parsable = True
            for prereq in value.prerequisites:
                if not self._parsable_quantities[prereq].has_files:
                    is_parsable = False
            value.is_parsable = is_parsable and value.has_files
            if quantity == value.nodeName:
                # This quantity and all of its alternatives represent a node.
                value.is_node = True
                for item in value.alternatives:
                    self._parsable_quantities[item].is_node = True
            # Add this quantity to the list of alternatives of another quantity.
            if value.is_alternative is not None:
                if value.is_alternative not in self._parsable_quantities:
                    # The quantity which this quantity is an alternative to is not in _parsable_quantities.
                    # Add an unparsable dummy quantity for it.
                    self.add_parsable_quantity(value.is_alternative, {})
                if quantity not in self._parsable_quantities[value.is_alternative].alternatives:
                    self._parsable_quantities[value.is_alternative].alternatives.append(quantity)

    def _check_and_validate_settings(self):
        """Check the settings and set which files should be parsed based on the input."""

        def add_quantity(quantity_to_add):
            """Check, whether a quantity or it's alternatives can be added."""
            for item in [quantity_to_add] + self._parsable_quantities[quantity_to_add].alternatives:
                if self._parsable_quantities[item].is_parsable:
                    self._quantities_to_parse.append(item)
                    return True
            return False

        for key, value in self._settings.iteritems():
            if not key.startswith('add_'):
                # only keys starting with 'add_' will change the behaviour of the parser so get the next one.
                continue
            if not value:
                # The quantity should not be added, so the corresponding files do not have to be parsed.
                continue
            quantity = key[4:]
            if quantity not in self._parsable_quantities:
                self.logger.warning('{0} has been requested by setting add_{0}'.format(quantity) +
                                    ' however it has not been implemented. Please check the docstrings' +
                                    ' in aiida_vasp.parsers.vasp.py for valid input.')
                continue

            # Found a node, which should be added, add it to the quantities to parse.
            # if all files required for this quantity have been retrieved. If there are
            # alternatives for this quantity also try those.
            success = add_quantity(quantity)

            if not success:
                # Neither the quantity nor it's alternatives could be added to the quantities_to_parse.
                # Gather a list of all the missing files and issue a warning.
                missing_files = []
                for quant in [quantity] + self._parsable_quantities[quantity].alternatives:
                    for missing_file in self._parsable_quantities[quant].missing_files:
                        missing_files.append(missing_file)

                missing_files = ", ".join(missing_files)
                self.logger.warning('{0} has been requested, however the following files'
                                    ''.format(quantity) +
                                    ' required for parsing have not been retrieved: {0}.'
                                    ''.format(missing_files))

    def _set_file_parsers(self):
        """
        Set the specific file parsers for OUTCAR, DOSCAR, EIGENVAL and vasprun.xml.

        Return False if a critical file is missing, which will abort the parsing.
        """

        for quantity in self._quantities_to_parse:
            for filename in self._parsable_quantities[quantity]['parsers']:
                if self._parsers[filename].parser is not None:
                    # This parser has already been checked.
                    continue
                file_to_parse = self.get_file(filename)
                self._parsers[filename].parser = self._parsers[filename]['parser_class'](
                    file_to_parse, calc_parser_cls=self)

    # pylint: disable=unused-argument, no-self-use
    @delegate()
    def get_quantity(self, quantity, settings):
        """Delegate to which the FileParsers will subscribe their get_quantity methods when they are initialised."""
        return {quantity: None}

    def get_inputs(self, quantity):
        """
        Return a quantity required as input for another quantity.

        This method will be called by the FileParsers in order to get a required input quantity
        from self._output_nodes. If the quantity is not in the dictionary the VaspParser will
        try to parse it. If a quantiy has been requested this way two times, parsing will be
        aborted because there is a cyclic dependency of the parsable items.
        """

        if quantity in self._requested_quantities:
            raise RuntimeError('{0} has been requested for parsing a second time.'.format(quantity) +
                               ' There is probably a cycle in the prerequisites of the' + ' parsable_items in the single FileParsers.')

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
