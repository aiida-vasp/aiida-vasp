"""
The standard stream parser interface for VASP.

----------------------------------------------
Contains the parsing interfaces to ``parsevasp`` used to parse standard streams
for VASP related notification, warnings and errors.
"""
# pylint: disable=abstract-method
import re

from parsevasp.stream import Stream
from aiida_vasp.parsers.content_parsers.base import BaseFileParser


class StreamParser(BaseFileParser):
    """Parser used for parsing errors and warnings from VASP."""

    OPEN_MODE = 'r'
    DEFAULT_SETTINGS = {'quantities_to_parse': ['notifications']}

    PARSABLE_QUANTITIES = {
        'notifications': {
            'inputs': [],
            'name': 'notifications',
            'prerequisites': [],
        }
    }

    def _init_from_handler(self, handler):
        """Initialize a ``parsevasp`` object of ``Stream`` using a file like handler.

        Parameters
        ----------
        handler : object
            A file like object that provides the necessary standard stream content to be parsed.

        """

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
        """Fetch the notifications that VASP generated.

        Returns
        -------
        notifications : list
            A list of all notifications from VASP. Each entry is a dict with the keys ``name``, ``kind``, ``message``
            and ``regex`` containing name of the message, what kind it is (``ERROR`` or ``WARNING``), a description
            of the notification and the regular expression detected as string values.

        """

        # ``parsevasp`` returns ``VaspStream`` objects, which we cannot serialize. We could serialize this, but
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
        """Fetch the errors that VASP generated.

        Returns
        -------
        errors : list
            A list of all errors from VASP. Each entry is a dict with the keys ``name``, ``kind``, ``message``
            and ``regex`` containing name of the message, what kind it is (``ERROR`` or ``WARNING``), a description
            of the error and the regular expression detected as string values.

        """

        return [item for item in self._content_parser.entries if item.kind == 'ERROR']

    @property
    def warnings(self):
        """Fetch the warnings that VASP generated.

        Returns
        -------
        warnings : list
            A list of all warnings from VASP. Each entry is a dict with the keys ``name``, ``kind``, ``message``
            and ``regex`` containing name of the message, what kind it is (``ERROR`` or ``WARNING``), a description
            of the error and the regular expression detected as string values.

        """

        return [item for item in self._content_parser.entries if item.kind == 'WARNING']

    @property
    def has_entries(self):
        """Check if there are notifications from VASP present according to the config after parsning.

        Returns
        -------
        entries : bool
            ``True`` if notifications was detected, ``False`` otherwise.

        """

        entries = self._content_parser.has_entries
        return entries

    @property
    def number_of_entries(self):
        """Find the number of unique notifications from VASP.

        Returns
        -------
        number_of_entries : int
            The number of unique notification entries that VASP generated.

        """

        number_of_entries = len(self._content_parser)
        return number_of_entries
