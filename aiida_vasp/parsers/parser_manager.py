"""Classes representing FileParsers in the VaspParser."""

from aiida_vasp.utils.extended_dicts import DictWithAttributes


class ParserManager(object):
    """
    A manager for FileParsers.

    :param vasp_parser: Reference to the VaspParser required for logging.
    :param quantities: Reference to a ParsableQuantities object for getting quantities.
    :param settings: A dictionary holding the 'parser_settings'.
    """

    def __init__(self, vasp_parser=None):
        self._parsers = {}
        self._quantities_to_parse = []

        self._vasp_parser = vasp_parser
        self._quantities = vasp_parser.quantities
        self._settings = vasp_parser.settings

        # Add all FileParsers from the requested set.
        for key, value in self._settings.parser_definitions.items():
            self.add_file_parser(key, value)

    def get_quantities_to_parse(self):
        return self._quantities_to_parse

    def remove(self, quantity):
        if quantity in self._quantities_to_parse:
            self._quantities_to_parse.remove(quantity)

    def add_file_parser(self, parser_name, parser_dict):
        """
        Add the definition of a fileParser to self._parsers.

        :param parser_name: Unique identifier of this parser. At the moment this coincides with the file name.
        :param parser_dict: Dict holding the FileParser definition.

        The required tags for the parser_dict can be found in PARSABLE_FILES. The FileParser must inherit
        from BaseFileParser and it will replace another previously defined fileParser with the same name.
        """

        parser_dict['parser'] = None
        parser_dict['quantities_to_parse'] = []
        self._parsers[parser_name] = DictWithAttributes(parser_dict)

    def add_quantity_to_parse(self, quantities):
        """Check, whether a quantity or it's alternatives can be added."""
        for quantity in quantities:
            if quantity.is_parsable:
                self._quantities_to_parse.append(quantity.name)
                return True
        return False

    def setup(self):

        self._set_quantities_to_parse()
        self._set_file_parsers()

    def _set_quantities_to_parse(self):
        """Set the quantities to parse list."""

        self._quantities_to_parse = []
        for quantity_name in self._settings.nodes:
            if not self._quantities.get_by_name(quantity_name):
                self._vasp_parser.logger.warning('{quantity} has been requested by setting '
                                                 'add_{quantity}. However it has not been implemented. '
                                                 'Please check the docstrings in aiida_vasp.parsers.vasp.py '
                                                 'for valid input.'.format(quantity=quantity_name))
                continue

            # Add this quantity or one of its alternatives to the quantities to parse.
            success = self.add_quantity_to_parse(self._quantities.get_equivalent_quantities(quantity_name))

            if not success:
                # Neither the quantity nor it's alternatives could be added to the quantities_to_parse.
                # Gather a list of all the missing files and issue a warning.
                missing_files = self._quantities.get_missing_files(quantity_name)
                self._vasp_parser.logger.warning('{quantity} has been requested, however the '
                                                 'following files required for parsing have not been '
                                                 'retrieved: {missing_files}.'.format(quantity=quantity_name, missing_files=missing_files))

    def _set_file_parsers(self):
        """Set the specific FileParsers."""

        for quantity in self._quantities_to_parse:
            for filename in self._quantities.get_by_name(quantity).parsers:
                if self._parsers[filename].parser is not None:
                    # This parser has already been checked, i.e. take the first
                    # available in the list that can be parsed (i.e. file exists)
                    continue
                file_to_parse = self._vasp_parser.get_file(filename)
                self._parsers[filename].parser = self._parsers[filename]['parser_class'](self._vasp_parser, file_path=file_to_parse)
