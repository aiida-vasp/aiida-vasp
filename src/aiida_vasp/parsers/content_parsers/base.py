"""
Base classes for the VASP content parsers.

===============================

Contains the base classes for the VASP content parsers.
"""
# pylint: disable=import-outside-toplevel

from __future__ import annotations

from typing import Any, BinaryIO, TextIO

from aiida.common import AIIDA_LOGGER
from aiida.orm import Data


class BaseFileParser:
    """
    Base class for all the content parsers which parse (read and write) VASP files.

    This class provides the interface for parsing and writing VASP files using parsevasp.
    It ensures integration with the AiiDA parsing framework and preparation before submission.
    Specific content parser interfaces should inherit from this class.

    There are two main usage patterns:

    1. **Parsing from file**: Initialize with ``handler`` (a file-like object), then fetch quantities
       using :meth:`get_quantity` with the appropriate key. The valid keys are defined in the
       ``PARSABLE_QUANTITIES`` class parameter of each subclass.

    2. **Writing from data**: Initialize with ``data`` (an AiiDA data node), then use :meth:`get_quantity`
       to retrieve the same node, or use :meth:`write` to write the content to a file.

    :param handler: File-like object containing content to be parsed. Used when parsing completed calculations.
    :type handler: file-like object, optional
    :param data: AiiDA data node. Used when writing VASP input files.
    :type data: object, optional
    :param settings: Parser settings, e.g. which quantities to compose into nodes.
    :type settings: dict, optional
    :param options: Parser options, e.g. extra options for content parsers such as selective dynamics for POSCAR.
    :type options: dict, optional
    """

    OPEN_MODE: str = 'r'
    PARSABLE_QUANTITIES: dict[str, Any] = {}
    DEFAULT_SETTINGS: dict[str, Any] = {'quantities_to_parse': []}

    def __init__(
        self,
        *,
        handler: TextIO | BinaryIO | None = None,
        data: Data | None = None,
        settings: dict[str, Any] | None = None,
        options: dict[str, Any] | None = None,
    ) -> None:
        super().__init__()
        # Make sure we only accept initialization with either ``handler`` or ``data``.
        if (handler is not None and data is not None) or (handler is None and data is None):
            raise TypeError('Supply at bare minimum either argument handler or data to initialize parser.')

        # Make sure logger messages in the parser are passed to the AiiDA logger.
        self._logger = AIIDA_LOGGER.getChild(self.__class__.__name__)

        # What quantities the specific content parser can provide.
        self._parsable_quantities = self.PARSABLE_QUANTITIES
        # The container for the parsed data when the ``get_quantity`` is executed, i.e. in the node composer
        # at a later stage.
        self._parsed_content = {}
        # The content parser, which will be an instance of one of the parsevasp parser classes.
        self._content_parser = None
        # Content data, which is an AiiDA data structure.
        self._content_data = None
        # Parser settings.
        self._set_settings(settings)
        # Parser options.
        self._options = options

        self.parser_notifications = {}

        # Set ``handler`` (parsing from some source) or ``data`` (eventually for example executing write)
        if handler is not None:
            self._init_from_handler(handler)
        if data is not None:
            # Check that the supplied handler is Data, one of the AiiDA supported data types.
            if isinstance(data, Data):
                self._init_from_data(data)
            else:
                raise TypeError('The supplied handler is not of Data type.')

    def get_all_quantities(self) -> tuple[dict[str, Any], dict[str, Any]]:
        """
        Fetch all quantities that can be parsed.

        :return: Tuple of (parsed, errored) dictionaries.
        :rtype: tuple
        """
        parsed = {}
        errored = {}
        for name in self.PARSABLE_QUANTITIES:
            try:
                parsed[name] = getattr(self, name)
            except (TypeError, ValueError, KeyError, AttributeError, IndexError) as error:
                errored[name] = error
                continue
        return parsed, errored

    @property
    def parsable_quantities(self) -> list[str]:
        """
        Fetch the quantities that this content parser can provide.

        :return: List of parsable quantity keys.
        :rtype: list
        """
        return self._parsable_quantities

    def _set_settings(self, settings: dict[str, Any] | None) -> None:
        """
        Set the settings to be used for the content parser.

        :param settings: The settings to be used for the content parser. Can be None if no settings
          is supplied to init. Defaults are then set.
        :type settings: None or dict
        """
        if settings is None:
            # Apply defaults
            self._settings = self.DEFAULT_SETTINGS
        else:
            self._settings = settings
        if not self._settings.get('quantities_to_parse'):
            # Did not find any quantities to parse in the settings, set it to the default
            # for each content parser or to an empty list of not defined
            self._settings['quantities_to_parse'] = (
                self.DEFAULT_SETTINGS['quantities_to_parse'] if self.DEFAULT_SETTINGS.get('quantities_to_parse') else []
            )
        # Let us make sure the quantities to parser is in a list form
        if not isinstance(self._settings.get('quantities_to_parse'), list):
            raise TypeError('The quantities_to_parse is not defined as a list of quantities.')

    def get_quantity(self, quantity_key: str) -> Any:
        """
        Fetch the required quantity from the content parser.

        Either fetch it from an existing AiiDA data structure, a parsed content dictionary
        if that exists, otherwise parse this specific quantity using the loaded instance,
        which is now a specific content parser.

        :param quantity_key: Key of the quantity to be fetched.
        :type quantity_key: str
        :return: The requested quantity, or None if not parsable.
        :rtype: object or None
        """
        if self._content_data is not None:
            # If we have already set an AiiDA data structure, return it.
            # This is straightforward in our case as there is for the PoscarParser,
            # KpointsParser a 1:1 mapping between the parser and the AiiDA data (if
            # we ignore conversions between representations etc.).
            return self._content_data
        # We continue assuming we need to parse this quantity
        if quantity_key not in self._parsable_quantities:
            # Check if this quantity can be parsed by this content parser, if not,
            # return None.
            return None
        if self._parsed_content.get(quantity_key) is None:
            # Parsed content does not contain this quantity,
            # most likely none of the content is parsed. Parse
            # all relevant content now and store.
            self._parsed_content = self._parse_content()

        return self._parsed_content.get(quantity_key)

    def write(self, path: str) -> None:
        """
        Write VASP content to file using the loaded content parser.

        :param path: Path to write the file to.
        :type path: str
        """

        if self._content_parser is None or self._content_data is None:
            # Only write if we have an AiiDA data structure or parser prepared.
            if self._content_parser is None:
                # If we do not have a parser loaded before write, we have an
                # AiiDA data structure. Make sure that is on the form parsevasp expects
                # by initializing the content parser instance.
                self._content_parser = self._content_data_to_content_parser()
            # Now a content parser should be ready and its content can be
            # written using parsevasp. But the content parser could still be None
            # if there is something with the data that can not be parsed.
            with open(path, 'w', encoding='utf8') as handler:
                self._content_parser.write(file_handler=handler)
        else:
            raise ValueError('The content parser has not been initialized or no AiiDA data structure is preparred.')

    def _init_from_handler(self, handler: TextIO | BinaryIO) -> None:
        """
        Initialize using a file-like object.

        Should be overridden in specific content parsers under ``content_parsers``
        if it will accept parsable content.

        :param handler: File-like object that provides the necessary content to be parsed.
        :type handler: object
        """

        raise NotImplementedError(f'{self.__class__.__name__} does not implement a _init_from_handler() method.')

    def _init_from_data(self, data: Data) -> None:
        """
        Initialize using an AiiDA data structure.

        Should be overridden in specific content parsers under ``content_parsers``
        if it will accept an AiiDA data structure. It should also check that the
        right structure is supplied.

        :param data: A valid AiiDA data structure object.
        :type data: object
        """

        raise NotImplementedError(f'{self.__class__.__name__} does not implement a _init_from_data() method.')

    def _content_data_to_content_parser(self) -> Any:
        """
        Convert an AiiDA data structure to a content parser instance relevant for that
        data structure. E.g. ``Poscar`` from ``parsevasp`` for an AiiDA ``StructureData``.

        Should be overridden in specific content parsers under ``content_parsers``
        if it will accept an AiiDA data structure. It should also check that the
        right structure is supplied.

        :return: Instance of a content parser from ``parsevasp``, e.g. ``Poscar``.
        :rtype: object
        """

        raise NotImplementedError(
            f'{self.__class__.__name__} does not implement a _content_data_to_content_parser() method.'
        )

    def _parse_content(self) -> dict[str, Any]:
        """
        Parse the quantities configured and parseable from the content.

        :return: Dictionary of parsed quantities.
        :rtype: dict
        """
        quantities_to_parse = self._settings.get('quantities_to_parse')
        result = {}
        if self._content_parser is None:
            # Parsevasp threw an exception, which means content could not be parsed.
            for quantity in quantities_to_parse:
                if quantity in self._parsable_quantities:
                    result[quantity] = None
            return result

        for quantity in quantities_to_parse:
            if quantity in self._parsable_quantities:
                # In case there is a - in the quantity, we assume we can
                # parse this quantity from multiple sources, remove source as we do not want to used
                # the source in the property name, i.e. use last element in the split
                quantity_splitted = quantity.split('-')
                quantity_splitted = quantity_splitted[-1]
                result[quantity] = getattr(self, quantity_splitted)

        return result
