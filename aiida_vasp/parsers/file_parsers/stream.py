"""
Standard stream parser for the VASP output.

--------------
A parser that analyses standard output for VASP related notification, warnings and
errors. Typically this is contained in the scheduler standard output.
"""

from parsevasp.stream import Stream
from aiida_vasp.parsers.file_parsers.parser import BaseFileParser, SingleFile
from aiida_vasp.parsers.error_handling import ErrorScanner

DEAULT_OPTIONS = {'quantities_to_parse': ['errors', 'warnings']}


class StreamParser(BaseFileParser):
    """Parser used for parsing errors and warnings from VASP."""

    PARSABLE_ITEMS = {
        'errors': {
            'inputs': [],
            'name': 'errors',
            'prerequisites': [],
        },
        'warnings': {
            'inputs': [],
            'name': 'warnings',
            'prerequisites': [],
        }
    }

    def __init__(self, *args, **kwargs):
        super(StreamParser, self).__init__(*args, **kwargs)
        self._notifications = None
        self.init_with_kwargs(**kwargs)

    def _init_with_file_path(self, path):
        """Init with a file path."""
        self._parsed_data = {}
        self.parsable_items = self.__class__.PARSABLE_ITEMS
        self._data_obj = SingleFile(path=path)

    def _init_with_data(self, data):
        """Init with SingleFileData."""
        self.parsable_items = self.__class__.PARSABLE_ITEMS
        self._init_with_file_path(data.get_file_abs_path())

    def _parse_file(self, inputs):
        """Parse the standard streams."""

        quantities_to_parse = DEAULT_OPTIONS.get('quantities_to_parse')
        if self.settings is not None:
            quantities_to_parse = self.settings.quantities_to_parse

        try:
            self._notifications = Stream(file_path=path, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abruptly while parsing the standard stream. Returning None.')
            self._notifications = None
        
        result = {}
        if notifications is None:
            # parsevasp threw an exception, which meas the standard stream could not be parsed.
            for quantity in quantities_to_parse:
                if quantity in self.parsable_items:
                    result[quantity] = None
            return result
        
        for quantity in quantities_to_parse:
            if quantity in self.parsable_items:
                result[quantity] = getattr(self, quantity)

        return result

    @property
    def errors(self):
        """Fetch the errors from parsevasp."""
        return self._notifications.entries['errors']

    @property
    def warnings(self):
        """Fetch the warnings from parsevasp."""
        return self._notifications.entries['warnings']
