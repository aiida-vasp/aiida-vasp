"""
Standard stream parser for the VASP output.

--------------
A parser that analyses standard output for VASP related notification, warnings and
errors. Typically this is contained in the scheduler standard output.
"""
import re

from parsevasp.stream import Stream
from aiida_vasp.parsers.file_parsers.parser import BaseFileParser, SingleFile

DEFAULT_OPTIONS = {'quantities_to_parse': ['notifications']}


class StreamParser(BaseFileParser):
    """Parser used for parsing errors and warnings from VASP."""

    PARSABLE_ITEMS = {
        'notifications': {
            'inputs': [],
            'name': 'notifications',
            'prerequisites': [],
        }
    }

    def __init__(self, *args, **kwargs):
        super(StreamParser, self).__init__(*args, **kwargs)
        self._stream = None
        self.init_with_kwargs(**kwargs)

    def _init_with_file_path(self, path):
        """Init with a file path."""
        self._parsed_data = {}
        self.parsable_items = self.__class__.PARSABLE_ITEMS
        self._data_obj = SingleFile(path=path)

        # Since the VASP output can be fairly large, we will parse it only
        # once and store the parsevasp Stream object.
        # First get any special config from the parser settings, else use the default
        stream_config = None
        history = False
        if self.settings is not None:
            stream_config = self.settings.get('stream_config', None)
            history = self.settings.get('stream_history', False)
        try:
            self._stream = Stream(file_path=path, logger=self._logger, history=history, config=stream_config)
        except SystemExit:
            self._logger.warning('Parsevasp exited abruptly when parsing the standard stream. Returning None.')
            self._stream = None

    def _init_with_data(self, data):
        """Init with SingleFileData."""
        self.parsable_items = self.__class__.PARSABLE_ITEMS
        self._init_with_file_path(data.get_file_abs_path())

    def _parse_file(self, inputs):
        """Parse the standard streams."""

        # Since all quantities will be returned by properties, we can't pass
        # inputs as a parameter, so we store them in self._parsed_data
        for key, value in inputs.items():
            self._parsed_data[key] = value

        quantities_to_parse = DEFAULT_OPTIONS.get('quantities_to_parse')
        if self.settings is not None and self.settings.quantities_to_parse:
            quantities_to_parse = self.settings.quantities_to_parse

        result = {}

        if self._stream is None:
            # parsevasp threw an exception, which means the standard stream could not be parsed.
            for quantity in quantities_to_parse:
                if quantity in self.parsable_items:
                    result[quantity] = None
            return result

        for quantity in quantities_to_parse:
            if quantity in self.parsable_items:
                result[quantity] = getattr(self, quantity)

        return result

    @property
    def notifications(self):
        """Fetch the notifications from parsevasp."""
        # Parsevasp returns VaspStream objects, which we cannot serialize. We could serialize this, but
        # eventually, we would like to move to a dedicated node for the notifications with its own data class.
        # This should be fixed in AiiDA core and coordinated across many plugins. For now, we convert the relevant info
        # into dict entries explicitly.
        notifications = []
        for item in self._stream.entries:
            if isinstance(item.regex, type(re.compile(''))):
                regex = item.regex.pattern
            else:
                regex = item.regex
            notifications.append({'name': item.shortname, 'kind': item.kind, 'message': item.message, 'regex': regex})

        return notifications

    @property
    def errors(self):
        """Fetch the errors from parsevasp."""
        return [item for item in self._stream.entries if item.kind == 'ERROR']

    @property
    def warnings(self):
        """Fetch the warnings from parsevasp."""
        return [item for item in self._stream.entries if item.kind == 'WARNING']

    @property
    def has_entries(self):
        """Fetch if there are streams present according to the config after parsning."""
        return self._stream.has_entries

    @property
    def number_of_entries(self):
        """Return a dict containing the number of unique streams detected."""
        return len(self._stream)
