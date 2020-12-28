"""
Parser manager.

---------------
Manages the parsing and executes the actual parsing by locating which file parsers
to use for each defined quantity.
"""


class ParserManager(object):  # pylint: disable=useless-object-inheritance
    """
    A manager for FileParsers.

    :param vasp_parser: Reference to the VaspParser required for logging.
    :param quantities: Reference to a ParsableQuantities object for getting quantities.
    :param settings: A dictionary holding the 'parser_settings'.
    """

    def __init__(self, vasp_parser_logger=None):
        self._vasp_parser_logger = vasp_parser_logger
        self._quantities_to_parse = []

    @property
    def quantities_to_parse(self):
        return self._quantities_to_parse

    def setup(self, quantities_to_parse=None, quantities=None, retrieve_list=None):
        """Collect quantity keys in quantities_to_parse"""

        self._quantities_to_parse = []
        for quantity_name in quantities_to_parse:
            if quantities.get_by_name(quantity_name):
                any_parsable = False
                for quantity in quantities.get_equivalent_quantities(quantity_name):
                    if quantity.is_parsable:
                        self._quantities_to_parse.append(quantity.original_name)
                        any_parsable = True
                if not any_parsable:
                    self._issue_warning(quantities.get_missing_files(quantity_name), retrieve_list, quantity_name)
            else:
                self._vasp_parser_logger.warning('{quantity} has been requested, '
                                                 'however its parser has not been implemented. '
                                                 'Please check the docstrings in aiida_vasp.parsers.vasp.py '
                                                 'for valid input.'.format(quantity=quantity_name))

    def _issue_warning(self, missing_files, retrieve_list, quantity_name):
        """
        Issue warning when no parsable quantity is found.

        Neither the quantity nor it's alternatives could be added to the quantities_to_parse.
        Gather a list of all the missing files and issue a warning.
        Check if the missing files are defined in the retrieve list

        """
        not_in_retrieve_list = None
        for item in missing_files:
            if item not in retrieve_list:
                not_in_retrieve_list = item
        self._vasp_parser_logger.warning('The quantity {quantity} has been requested for parsing, however the '
                                         'following files required for parsing it have not been '
                                         'retrieved: {missing_files}.'.format(quantity=quantity_name, missing_files=missing_files))
        if not_in_retrieve_list is not None:
            self._vasp_parser_logger.warning(
                'The file {not_in_retrieve_list} is not present '
                'in the list of files to be retrieved. If you want to add additional '
                'files, please make sure to define it in the ADDITIONAL_RETRIEVE_LIST, '
                'which is an option given to calculation settings.'.format(not_in_retrieve_list=not_in_retrieve_list))
