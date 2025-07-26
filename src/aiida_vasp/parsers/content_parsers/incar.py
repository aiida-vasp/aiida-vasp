"""
The ``INCAR`` parser interface.

Contains the parsing interfaces to parsevasp used to parse ``INCAR`` content.
"""

from aiida import orm
from parsevasp.incar import Incar

from aiida_vasp.parsers.content_parsers.base import BaseFileParser


class IncarParser(BaseFileParser):
    """The parser interface that enables parsing of ``INCAR`` content.

    The parser is triggered by using the ``incar`` quantity key.
    """

    DEFAULT_SETTINGS = {'quantities_to_parse': ['incar']}

    PARSABLE_QUANTITIES = {
        'incar': {'inputs': [], 'name': 'incar', 'prerequisites': []},
    }

    def __init__(self, *args, validate_tags=True, **kwargs):
        self._validate_tags = validate_tags
        super().__init__(*args, **kwargs)

    def _init_from_handler(self, handler):
        """Initialize a ``parsevasp`` object of ``Incar`` using a file like handler.

        :param handler: A file like object that provides the necessary ``INCAR`` content to be parsed
        :type handler: object
        """

        try:
            self._content_parser = Incar(file_handler=handler, logger=self._logger, validate_tags=self._validate_tags)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    def _init_from_data(self, data):
        """Initialize using an AiiDA ``Dict`` instance.

        :param data: A valid AiiDA ``Dict`` object
        :type data: object
        :raises TypeError: If the supplied AiiDA data structure is not a Dict
        """

        if isinstance(data, orm.Dict):
            self._content_data = data
        else:
            raise TypeError('The supplied AiiDA data structure is not a Dict.')

    @property
    def incar(self):
        """Return the parameters in the ``INCAR``.

        :returns: A dictionary containing the parameter tags as keys and its settings as values.
            ``None`` is returned if the quantity can not be parsed.
        :rtype: dict or None
        """
        if self._content_parser is not None:
            params = self._content_parser.get_dict()
            return params
        return None

    def _content_data_to_content_parser(self):
        """Convert an AiiDA ``Dict`` to a content parser instance of ``Incar`` from ``parsevasp``.

        :returns: An instance of ``Incar`` from ``parsevasp``
        :rtype: object
        """
        # Filter away None values from the dictionary - these are not valid for ``parsevasp``
        # This allow easier workflow control and parameters merging - setting a key to None means it should not be
        # set
        dictionary = {key: value for key, value in self._content_data.get_dict().items() if value is not None}

        # We brake hard if ``parsevasp`` fail here. If we can not write we will not try another parser.
        content_parser = Incar(incar_dict=dictionary, logger=self._logger, validate_tags=self._validate_tags)

        return content_parser
