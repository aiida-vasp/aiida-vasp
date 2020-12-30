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
    self.name corresponds to link name.

    """

    def __init__(self, original_name, init_dict):
        self.original_name = original_name
        self.alternatives = []
        self.prerequisites = []
        self.inputs = []
        self.name = None
        self.file_name = None
        super(ParsableQuantity, self).__init__(init_dict)

        if self.get('name') is None:
            self.name = original_name


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
        self._missing_filenames = {}
        self._quantity_keys = None
        self._equiv_node_keys = {}
        self._waiting_quantity_items = {}
        self._quantity_items = {}
        self._quantity_keys_to_filenames = {}

    @property
    def quantity_keys_to_parse(self):
        return self._quantity_keys_to_parse

    def add_parsable_quantity(self, quantity_key, quantity_dict):
        self._waiting_quantity_items[quantity_key] = quantity_dict

    def remove_parsable_quantity(self, quantity_key):
        _ = self._parsable_quantities.pop(quantity_key, None)

    def get_equivalent_quantities(self, quantity_name):
        """Get a list of equivalent quantities."""
        return [self._parsable_quantities[alt_quantity_key] for alt_quantity_key in self._alternatives[quantity_name]]

    def get_by_name(self, quantity_key):
        """Get a quantity by name."""
        return self._parsable_quantities.get(quantity_key)

    def setup(self, retrieved_filenames=None, parser_definitions=None, quantity_keys_to_parse=None):
        """Set the parsable_quantities dictionary based on parsable_items obtained from the FileParsers."""
        self._add_quantities_in_definitions(parser_definitions)
        self._create_containers_of_equiv_node_keys()
        self._show(self._quantity_keys)
        self._put_alternatives_in_equiv_node_keys()
        self._show(self._equiv_node_keys)
        self._set_missing_filenames(retrieved_filenames)
        self._debug_show(self._missing_filenames)
        self._check_parsable(retrieved_filenames)
        self._screen_quantity_keys_to_parse(quantity_keys_to_parse, retrieved_filenames)

    def _show(self, x):
        print('---debug---')
        if isinstance(x, dict):
            for key, value in x.items():
                print(key, value)
        elif isinstance(x, list):
            for value in x:
                print(value)
        print('---debug---')

    def _debug_show(self, x):
        self._show(x)
        raise

    def _add_quantities_in_definitions(self, parser_definitions):
        """Gather all quantity keys in definitions of the file parsers."""
        for filename, value in parser_definitions.items():
            for quantity_key, quantity_dict in value['parser_class'].PARSABLE_ITEMS.items():
                self._add_quantity(quantity_key, quantity_dict, filename)

        for quantity_key, quantity_dict in self._waiting_quantity_items.items():
            if 'file_name' in quantity_dict:
                self._add_quantity(quantity_key, quantity_dict, quantity_dict['file_name'])
            else:
                raise RuntimeError('The added quantity {quantity} has to contain file_name entry.'.format(quantity=quantity_key))

        self._quantity_keys = list(self._quantity_items.keys())

    def _add_quantity(self, quantity_key, quantity_dict, filename):
        if quantity_key in self._quantity_items:
            # This quantity has already been added so it is not unique.
            raise RuntimeError('The quantity {quantity} defined in {filename} has been '
                               'defined by two FileParser classes. Quantity names must '
                               'be unique. If both quantities are equivalent, define one '
                               'as an alternative for the other.'.format(quantity=quantity_key, filename=filename))
        self._quantity_keys_to_filenames[quantity_key] = filename
        self._quantity_items[quantity_key] = quantity_dict

    def _create_containers_of_equiv_node_keys(self):
        """
        Create containars of quantity keys to equivalent node keys

        Put the first entry for each node key.
        The first entry of each container is that with quantity_key == quantity['name'] if found.

        """
        for quantity_key in self._quantity_keys:
            quantity_dict = self._quantity_items[quantity_key]
            if 'name' in quantity_dict:
                name = quantity_dict['name']
            else:
                name = quantity_key
            if name not in self._equiv_node_keys:
                self._equiv_node_keys[name] = []
            if quantity_key == name and name not in self._equiv_node_keys[name]:
                self._equiv_node_keys[name].append(quantity_key)

        for quantity_key in self._quantity_keys:
            quantity_dict = self._quantity_items[quantity_key]
            if 'name' in quantity_dict:
                name = quantity_dict['name']
                if quantity_key not in self._equiv_node_keys[name]:
                    self._equiv_node_keys[name].append(quantity_key)

    def _put_alternatives_in_equiv_node_keys(self):
        """Put alternative node keys to equiv_node_keys"""

        for quantity_key in self._quantity_keys:
            quantity_dict = self._quantity_items[quantity_key]
            for node_key in self._equiv_node_keys:
                if 'name' in quantity_dict:
                    name = quantity_dict['name']
                else:
                    name = quantity_key
                if name == node_key:
                    if 'alternatives' in quantity_dict:
                        for alt_node_key in quantity_dict['alternatives']:
                            if alt_node_key not in self._equiv_node_keys[node_key]:
                                self._equiv_node_keys[node_key].append(alt_node_key)

    def _set_missing_filenames(self, retrieved_filenames):
        """Collect lists of all missing files for quantities."""
        for quantity_key in self._quantity_keys:
            missing_filenames = []
            filename = self._quantity_keys_to_filenames[quantity_key]
            if filename not in retrieved_filenames:
                if filename not in missing_filenames:
                    missing_filenames.append(filename)
            self._missing_filenames[quantity_key] = missing_filenames

    def _check_parsable(self, retrieved_filenames):
        """Check for every quantity, whether all of the prerequisites are parsable."""
        for quantity_key in self._parsable_quantities:
            quantity = self._parsable_quantities[quantity_key]
            is_parsable = True
            for prereq in quantity.prerequisites:
                if self._parsable_quantities.get(prereq) is None:
                    is_parsable = False
                elif self._missing_filenames[prereq]:
                    is_parsable = False
            quantity.is_parsable = is_parsable and not self._missing_filenames[quantity_key]

    def _screen_quantity_keys_to_parse(self, quantity_keys_to_parse, retrieve_filenames):
        """Collect quantity keys in quantity_keys_to_parse"""
        self._quantity_keys_to_parse = []
        for quantity_key in quantity_keys_to_parse:
            if self.get_by_name(quantity_key):
                any_parsable = False
                quantity_name = self._parsable_quantities[quantity_key].name
                for quantity in self.get_equivalent_quantities(quantity_name):
                    if quantity.is_parsable:
                        self._quantity_keys_to_parse.append(quantity_key)
                        any_parsable = True
                        break
                if not any_parsable:
                    self._issue_warning(retrieve_filenames, quantity_key)
            else:
                self._vasp_parser_logger.warning('{quantity} has been requested, '
                                                 'however its parser has not been implemented. '
                                                 'Please check the docstrings in aiida_vasp.parsers.vasp.py '
                                                 'for valid input.'.format(quantity=quantity_key))

    def _issue_warning(self, retrieve_filenames, quantity_key):
        """
        Issue warning when no parsable quantity is found.

        Neither the quantity nor it's alternatives could be added to the quantity_keys_to_parse.
        Gather a list of all the missing files and issue a warning.
        Check if the missing files are defined in the retrieve list

        """
        missing_files = self._missing_filenames[quantity_key]
        not_in_retrieve_filenames = None
        for item in missing_files:
            if item not in retrieve_filenames:
                not_in_retrieve_filenames = item
        self._vasp_parser_logger.warning('The quantity {quantity} has been requested for parsing, however the '
                                         'following files required for parsing it have not been '
                                         'retrieved: {missing_files}.'.format(quantity=quantity_key, missing_files=missing_files))
        if not_in_retrieve_filenames is not None:
            self._vasp_parser_logger.warning(
                'The file {not_in_retrieve_filenames} is not present '
                'in the list of retrieved files. If you want to add additional '
                'files, please make sure to define it in the ADDITIONAL_RETRIEVE_LIST, '
                'which is an option given to calculation settings.'.format(not_in_retrieve_filenames=not_in_retrieve_filenames))
