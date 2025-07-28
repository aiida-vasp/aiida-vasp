"""
The ``CHGCAR`` parser interface.

Contains the parsing interfaces to ``parsevasp`` used to parse ``CHGCAR`` content.
"""

# pylint: disable=abstract-method
from typing import TextIO

from parsevasp.chgcar import Chgcar

from aiida_vasp.parsers.content_parsers.base import BaseFileParser


class ChgcarParser(BaseFileParser):
    """The parser interface that enables parsing of ``CHGCAR`` content.

    The parser is triggered by using the ``charge_density`` and/or ``magnetization_density`` quantity key.

    """

    DEFAULT_SETTINGS = {'quantities_to_parse': ['charge_density']}

    PARSABLE_QUANTITIES = {
        'charge_density': {'inputs': [], 'name': 'charge_density', 'prerequisites': []},
        'magnetization_density': {'inputs': [], 'name': 'magnetization_density', 'prerequisites': []},
    }

    def _init_from_handler(self, handler: TextIO) -> None:
        """Initialize a ``parsevasp`` object of ``Chgcar`` using a file like handler.

        :param handler: A file like object that provides the necessary ``CHGCAR`` content to be parsed.
        :type handler: file-like object
        """

        try:
            self._content_parser = Chgcar(file_handler=handler, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    @property
    def charge_density(self):
        """Return the charge density.

        :returns: A NumPy array containing the charge density in the unit cell in C order.
        :rtype: ndarray
        """

        return self._content_parser.charge_density

    @property
    def magnetization_density(self):
        """Return the magnetization density.

        :returns: If collinear spin calculations have been performed, a NumPy array containing
                  the magnetization density in the unit cell in C order is returned. If however
                  a non-collinear spin calculation have been performed, a dictionary is returned
                  with keys `x`, `y` and `z`, each containing the same NumPy array, but for each
                  direction.
        :rtype: dict or ndarray
        """
        return self._content_parser.magnetization_density
