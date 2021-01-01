"""
Parser quantity configuration.

------------------------------
Contains the representation of quantities that users want to parse.
"""


class ParsableQuantities(object):  # pylint: disable=useless-object-inheritance
    """
    A Database of parsable quantities.

    - Get the parsable quantities from the FileParsers and initialise them.
    - Provide ways to get parsable quantities.
    """

    def __init__(self, vasp_parser_logger=None):
        self._vasp_parser_logger = vasp_parser_logger

        self._quantity_keys_to_parse = None
        self._equiv_quantity_keys = None
        self._quantity_keys_to_filenames = None

        self._missing_filenames = None
        self._quantity_items = None
        self._waiting_quantity_items = {}

    @property
    def quantity_keys_to_parse(self):
        """List of quantity keys to parse after screening"""
        return self._quantity_keys_to_parse

    @property
    def quantity_keys_to_filenames(self):
        """Dictionary of quantity key -> file name"""
        return self._quantity_keys_to_filenames

    @property
    def equivalent_quantity_keys(self):
        """
        Equivalent quantity keys of quantity names

        Returns
        -------
        dict
            {quantity_name: [quantity_key, ...], ...}

        """
        return self._equiv_quantity_keys

    def add_parsable_quantity(self, quantity_key, quantity_dict):
        """Put parsable quantity in the waiting list"""
        self._waiting_quantity_items[quantity_key] = quantity_dict

    def setup(self, retrieved_filenames=None, parser_definitions=None, quantity_names_to_parse=None):
        """Set the parsable_quantities dictionary based on parsable_items obtained from the FileParsers."""

        def _show(var, var_name):
            print('---%s ---' % var_name)
            if isinstance(var, dict):
                for key, value in var.items():
                    print(key, value)
            elif isinstance(var, list):
                for value in var:
                    print(value)
            else:
                print(var)
            print('---%s ---' % var_name)

        show_screening_steps = False

        self._quantity_items, self._quantity_keys_to_filenames = self._get_quantity_items_from_definitions(parser_definitions)
        if show_screening_steps:
            _show(self._quantity_items, 'self._quantity_items')
            _show(self._quantity_keys_to_filenames, 'self._quantity_keys_to_filenames')

        self._equiv_quantity_keys = self._create_containers_of_equiv_quantity_keys()
        if show_screening_steps:
            _show(self._equiv_quantity_keys, 'self._equiv_quantity_keys')

        self._missing_filenames = self._identify_missing_filenames(retrieved_filenames, parser_definitions.keys())
        if show_screening_steps:
            _show(retrieved_filenames, 'retrieved_filenames')
            _show(parser_definitions.keys(), 'parser_definitions.keys()')
            _show(self._missing_filenames, 'self._missing_filenames')

        parsable_quantity_keys = self._get_parsable_quantity_keys()
        if show_screening_steps:
            _show(parsable_quantity_keys, 'parsable_quantity_keys')

        self._quantity_keys_to_parse = self._get_quantity_keys_to_parse(parsable_quantity_keys, quantity_names_to_parse,
                                                                        retrieved_filenames)
        if show_screening_steps:
            _show(quantity_names_to_parse, 'quantity_names_to_parse')
            _show(self._quantity_keys_to_parse, 'self._quantity_keys_to_parse')

    def _get_quantity_items_from_definitions(self, parser_definitions):
        """
        Gather all quantity keys in definitions of the file parsers.

        quantity_key has to be unique.

        """
        _quantity_keys_to_filenames = {}
        _quantity_items = {}

        for filename, value in parser_definitions.items():
            for quantity_key, quantity_dict in value['parser_class'].PARSABLE_ITEMS.items():
                self._add_quantity(_quantity_items, _quantity_keys_to_filenames, quantity_key, quantity_dict, filename)

        for quantity_key, quantity_dict in self._waiting_quantity_items.items():
            if 'file_name' in quantity_dict:
                self._add_quantity(_quantity_items, _quantity_keys_to_filenames, quantity_key, quantity_dict, quantity_dict['file_name'])
            else:
                raise RuntimeError('The added quantity {quantity} has to contain file_name entry.'.format(quantity=quantity_key))

        return _quantity_items, _quantity_keys_to_filenames

    def _add_quantity(self, _quantity_items, _quantity_keys_to_filenames, quantity_key, quantity_dict, filename):  # pylint: disable=no-self-use
        """
        Helper to store quantity_dict in self._quantity_items

        Note
        ----
        When 'name' key is unavailable in quantity_dict,
        quantity_dict['name'] = quantity_key is set.

        """
        if quantity_key in _quantity_items:
            raise RuntimeError('The quantity {quantity} defined in {filename} has been '
                               'defined by two FileParser classes. Quantity names must '
                               'be unique. If both quantities are equivalent, define one '
                               'as an alternative for the other.'.format(quantity=quantity_key, filename=filename))
        _quantity_keys_to_filenames[quantity_key] = filename
        _quantity_dict = quantity_dict.copy()
        if 'name' not in _quantity_dict:
            _quantity_dict['name'] = quantity_key
        _quantity_items[quantity_key] = _quantity_dict

    def _create_containers_of_equiv_quantity_keys(self):
        """Create containars of quantity keys equivalent to the same quantity names"""

        _equiv_quantity_keys = {}

        for quantity_key, quantity_dict in self._quantity_items.items():
            name = quantity_dict['name']
            if name not in _equiv_quantity_keys:
                _equiv_quantity_keys[name] = []
            if quantity_key == name and name not in _equiv_quantity_keys[name]:
                _equiv_quantity_keys[name].append(quantity_key)

        for quantity_key, quantity_dict in self._quantity_items.items():
            quantity_dict = self._quantity_items[quantity_key]
            name = quantity_dict['name']
            if quantity_key not in _equiv_quantity_keys[name]:
                _equiv_quantity_keys[name].append(quantity_key)

        for quantity_key, quantity_dict in self._quantity_items.items():
            if 'alternatives' in quantity_dict:
                name = quantity_dict['name']
                if name in _equiv_quantity_keys:
                    for alt_quantity_key in quantity_dict['alternatives']:
                        if alt_quantity_key not in _equiv_quantity_keys[name]:
                            _equiv_quantity_keys[name].append(alt_quantity_key)

        return _equiv_quantity_keys

    def _identify_missing_filenames(self, retrieved_filenames, filenames_in_parser_definitions):
        """Identify missing files for quantities."""
        _missing_filenames = {}
        for quantity_key in self._quantity_items:
            filename = self._quantity_keys_to_filenames[quantity_key]
            if filename not in retrieved_filenames or filename not in filenames_in_parser_definitions:
                _missing_filenames[quantity_key] = filename
        return _missing_filenames

    def _get_parsable_quantity_keys(self):
        """Check for every quantity, whether all of the prerequisites are parsable."""
        _parsable_quantity_keys = []
        for quantity_key, quantity_dict in self._quantity_items.items():
            is_parsable = True
            if 'prerequisites' in quantity_dict:
                for prereq in quantity_dict['prerequisites']:
                    if prereq not in self._quantity_items or prereq in self._missing_filenames:
                        is_parsable = False
            if is_parsable and quantity_key not in self._missing_filenames:
                _parsable_quantity_keys.append(quantity_key)
        return _parsable_quantity_keys

    def _get_quantity_keys_to_parse(self, parsable_quantity_keys, quantity_names_to_parse, retrieve_filenames):
        """Collect quantity_names_to_parse"""
        _quantity_keys_to_parse = []
        for quantity_name in quantity_names_to_parse:
            if quantity_name in self._equiv_quantity_keys:
                is_parsable = False
                for quantity_key in self._equiv_quantity_keys[quantity_name]:
                    if quantity_key in parsable_quantity_keys:
                        _quantity_keys_to_parse.append(quantity_key)
                        is_parsable = True
                if not is_parsable:
                    self._issue_warning(retrieve_filenames, quantity_name)
            else:
                self._vasp_parser_logger.warning('{quantity} has been requested, '
                                                 'however its parser has not been implemented. '
                                                 'Please check the docstrings in aiida_vasp.parsers.vasp.py '
                                                 'for valid input.'.format(quantity=quantity_name))
        return _quantity_keys_to_parse

    def _issue_warning(self, retrieve_filenames, quantity_name):
        """
        Issue warning when no parsable quantity is found.

        Neither the quantity nor it's alternatives could be added to the quantity_keys_to_parse.
        Gather a list of all the missing files and issue a warning.
        Check if the missing files are defined in the retrieve list

        """
        missing_files = self._missing_filenames[quantity_name]
        not_in_retrieve_filenames = None
        for item in missing_files:
            if item not in retrieve_filenames:
                not_in_retrieve_filenames = item
        self._vasp_parser_logger.warning('The quantity {quantity} has been requested for parsing, however the '
                                         'following files required for parsing it have not been '
                                         'retrieved: {missing_files}.'.format(quantity=quantity_name, missing_files=missing_files))
        if not_in_retrieve_filenames is not None:
            self._vasp_parser_logger.warning(
                'The file {not_in_retrieve_filenames} is not present '
                'in the list of retrieved files. If you want to add additional '
                'files, please make sure to define it in the ADDITIONAL_RETRIEVE_LIST, '
                'which is an option given to calculation settings.'.format(not_in_retrieve_filenames=not_in_retrieve_filenames))
