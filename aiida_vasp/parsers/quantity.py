"""
Parser quantity configuration.

------------------------------
Contains the representation of quantities that users want to parse.
"""
# pylint: disable=import-outside-toplevel
from aiida_vasp.utils.extended_dicts import DictWithAttributes


class ParsableQuantity(DictWithAttributes):
    """
    Container class for parsable quantities.

    Example from PARSABLE_ITEMS in VasprunParser:
        self._parsable_items = {'structure': {
            'inputs': [],
            'name': 'structure',
            'prerequisites': [],
            'alternatives': ['poscar-structure']
        } ... }

        key = 'structure'
        value = items['structure']
        parsable_quantity = ParsableQuantity(key, value, None)

    self.file_name may be originated from keys of FILE_PARSER_SETS
    files_list may come from parser._retrieved_content.keys(), i.e.,
    [retrieved_file.name for retrieved_file in parser.retrieved.list_objects()].

    """

    def __init__(self, name, init_dict, files_list):
        self.original_name = name

        # The following items are given from init_dict.
        self.alternatives = []
        self.prerequisites = []
        self.inputs = []
        self.name = None
        super(ParsableQuantity, self).__init__(init_dict)

        if self.get('name') is None:
            self.name = name
        if self.get('file_name') is None:
            self.file_name = 'MISSING_FILE_NAME'

        # Check whether the file required for parsing this quantity have been retrieved.
        self.missing_files = []
        if files_list is None or self.file_name not in files_list:
            self.missing_files.append(self.file_name)


class ParsableQuantities(object):  # pylint: disable=useless-object-inheritance
    """
    A Database of parsable quantities.

    - Get the parsable quantities from the FileParsers and initialise them.
    - Provide ways to get parsable quantities.
    """

    def __init__(self, vasp_parser_logger=None):
        self._parsable_quantities = {}
        self._vasp_parser_logger = vasp_parser_logger
        self._quantity_keys_to_parse = []

    @property
    def quantity_keys_to_parse(self):
        return self._quantity_keys_to_parse

    def add_parsable_quantity(self, quantity_key, quantity_dict, retrieved_files=None):
        self._parsable_quantities[quantity_key] = ParsableQuantity(quantity_key, quantity_dict, retrieved_files)

    def remove_parsable_quantity(self, quantity_key):
        _ = self._parsable_quantities.pop(quantity_key, None)

    def get_equivalent_quantities(self, quantity_key):
        """Get a list of equivalent quantities."""
        quantity = self._parsable_quantities[quantity_key]
        return [quantity] + [self._parsable_quantities[alternative_quantity_key] for alternative_quantity_key in quantity.alternatives]

    def get_by_name(self, quantity_key):
        """Get a quantity by name."""
        return self._parsable_quantities.get(quantity_key)

    def get_missing_files(self, quantity_key):
        """Return a list with all missing files for a quantity."""
        missing_files = []
        for quantity in self.get_equivalent_quantities(quantity_key):
            for missing_file in quantity.missing_files:
                missing_files.append(missing_file)

        return missing_files

    def setup(self, retrieved_filenames=None, parser_definitions=None):
        """Set the parsable_quantities dictionary based on parsable_items obtained from the FileParsers."""

        self._add_parsable(retrieved_filenames, parser_definitions)
        self._collect_alternatives()

    def _add_parsable(self, retrieved_filenames, parser_definitions):
        """Gather all parsable items as defined in the file parsers."""

        self._parsable_quantities = {}
        for file_name, value in parser_definitions.items():
            for quantity_key, quantity_dict in value['parser_class'].PARSABLE_ITEMS.items():
                if quantity_key in self._parsable_quantities:
                    # This quantity has already been added so it is not unique.
                    raise RuntimeError('The quantity {quantity} defined in {filename} has been '
                                       'defined by two FileParser classes. Quantity names must '
                                       'be unique. If both quantities are equivalent, define one '
                                       'as an alternative for the other.'.format(quantity=quantity_key, filename=file_name))
                # Create quantity objects.
                quantity_dict['file_name'] = file_name
                self.add_parsable_quantity(quantity_key, quantity_dict, retrieved_filenames)

    def _collect_alternatives(self):
        """Check the consistency and collect alternatives."""

        parsable_quantity_keys = list(self._parsable_quantities.keys())

        # Setup all alternatives:
        for quantity_key in parsable_quantity_keys:
            quantity = self._parsable_quantities[quantity_key]
            if quantity_key != quantity.name:
                # This quantity is considered as an alternative to quantity.name.
                if quantity.name not in self._parsable_quantities:
                    # Add a dummy quantity for it.
                    self.add_parsable_quantity(quantity.name, {})
                if quantity_key not in self._parsable_quantities[quantity.name].alternatives:
                    self._parsable_quantities[quantity.name].alternatives.append(quantity_key)

        # Check for every quantity, whether all of the prerequisites are parsable.
        for quantity_key in self._parsable_quantities:
            quantity = self._parsable_quantities[quantity_key]
            is_parsable = True

            for prereq in quantity.prerequisites:
                if self._parsable_quantities.get(prereq) is None:
                    is_parsable = False
                    continue
                if self._parsable_quantities[prereq].missing_files:
                    is_parsable = False

            quantity.is_parsable = is_parsable and not quantity.missing_files

    def screen_quantity_keys_to_parse(self, quantity_keys_to_parse=None, retrieve_list=None):
        """Collect quantity keys in quantity_keys_to_parse"""

        self._quantity_keys_to_parse = []
        for quantity_key in quantity_keys_to_parse:
            if self.get_by_name(quantity_key):
                any_parsable = False
                for quantity in self.get_equivalent_quantities(quantity_key):
                    if quantity.is_parsable:
                        self._quantity_keys_to_parse.append(quantity.original_name)
                        any_parsable = True
                if not any_parsable:
                    self._issue_warning(self.get_missing_files(quantity_key), retrieve_list, quantity_key)
            else:
                self._vasp_parser_logger.warning('{quantity} has been requested, '
                                                 'however its parser has not been implemented. '
                                                 'Please check the docstrings in aiida_vasp.parsers.vasp.py '
                                                 'for valid input.'.format(quantity=quantity_key))

    def _issue_warning(self, missing_files, retrieve_list, quantity_key):
        """
        Issue warning when no parsable quantity is found.

        Neither the quantity nor it's alternatives could be added to the quantity_keys_to_parse.
        Gather a list of all the missing files and issue a warning.
        Check if the missing files are defined in the retrieve list

        """
        not_in_retrieve_list = None
        for item in missing_files:
            if item not in retrieve_list:
                not_in_retrieve_list = item
        self._vasp_parser_logger.warning('The quantity {quantity} has been requested for parsing, however the '
                                         'following files required for parsing it have not been '
                                         'retrieved: {missing_files}.'.format(quantity=quantity_key, missing_files=missing_files))
        if not_in_retrieve_list is not None:
            self._vasp_parser_logger.warning(
                'The file {not_in_retrieve_list} is not present '
                'in the list of files to be retrieved. If you want to add additional '
                'files, please make sure to define it in the ADDITIONAL_RETRIEVE_LIST, '
                'which is an option given to calculation settings.'.format(not_in_retrieve_list=not_in_retrieve_list))
