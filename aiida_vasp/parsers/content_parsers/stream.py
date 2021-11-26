"""
Standard stream parser for the VASP output.

--------------
Contains the parsing interfaces to parsevasp used to parse standard streams
for VASP related notification, warnings and errors.
"""
import re

from parsevasp.stream import Stream
from aiida_vasp.parsers.content_parsers.base import BaseFileParser


class StreamParser(BaseFileParser):
    """Parser used for parsing errors and warnings from VASP."""

    DEFAULT_OPTIONS = {'quantities_to_parse': ['notifications']}

    PARSABLE_QUANTITIES = {
        'notifications': {
            'inputs': [],
            'name': 'notifications',
            'prerequisites': [],
        }
    }

    def _init_from_handler(self, handler):
        """Initialize using a file like handler."""
        # First get any special config from the parser settings, else use the default
        stream_config = None
        history = False
        if self._settings is not None:
            stream_config = self._settings.get('stream_config', None)
            history = self._settings.get('stream_history', False)

        try:
            self._content_parser = Stream(file_handler=handler, logger=self._logger, history=history, config=stream_config)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    @property
    def notifications(self):
        """Fetch the notifications from parsevasp."""
        # Parsevasp returns VaspStream objects, which we cannot serialize. We could serialize this, but
        # eventually, we would like to move to a dedicated node for the notifications with its own data class.
        # This should be fixed in AiiDA core and coordinated across many plugins. For now, we convert the relevant info
        # into dict entries explicitly.
        notifications = []
        for item in self._content_parser.entries:
            if isinstance(item.regex, type(re.compile(''))):
                regex = item.regex.pattern
            else:
                regex = item.regex
            notifications.append({'name': item.shortname, 'kind': item.kind, 'message': item.message, 'regex': regex})

        return notifications

    @property
    def errors(self):
        """Fetch the errors from parsevasp."""
        return [item for item in self._content_parser.entries if item.kind == 'ERROR']

    @property
    def warnings(self):
        """Fetch the warnings from parsevasp."""
        return [item for item in self._content_parser.entries if item.kind == 'WARNING']

    @property
    def has_entries(self):
        """Fetch if there are streams present according to the config after parsning."""
        return self._content_parser.has_entries

    @property
    def number_of_entries(self):
        """Return a dict containing the number of unique streams detected."""
        return len(self._content_parser)
