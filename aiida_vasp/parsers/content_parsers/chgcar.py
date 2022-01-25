"""
The ``CHGCAR`` parser interface.

--------------------------------
Contains the parsing interfaces to ``parsevasp`` used to parse ``CHGCAR`` content.
"""
# pylint: disable=abstract-method
from aiida_vasp.parsers.content_parsers.base import BaseFileParser

from parsevasp.chgcar import Chgcar


class ChgcarParser(BaseFileParser):
    """The parser interface that enables parsing of ``CHGCAR`` content.

    The parser is triggered by using the ``charge_density`` and/or ``magnetization_density`` quantity key.

    """

    DEFAULT_SETTINGS = {'quantities_to_parse': ['charge_density']}

    PARSABLE_QUANTITIES = {
        'charge_density': {
            'inputs': [],
            'name': 'charge_density',
            'prerequisites': []
        },
        'magnetization_density': {
            'inputs': [],
            'name': 'magnetization_density',
            'prerequisites': []
        }
    }

    def _init_from_handler(self, handler):
        """Initialize a ``parsevasp`` object of ``Chgcar`` using a file like handler.

        Parameters
        ----------
        handler : object
            A file like object that provides the necessary ``CHGCAR`` content to be parsed.

        """

        try:
            self._content_parser = Chgcar(file_handler=handler, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    @property
    def charge_density(self):
        """
        Return the charge density.

        Returns
        -------
        charge_density : ndarray
            A NumPy array containing the charge density in the unit cell in C order.

        """

        charge_density = self._content_parser.charge_density
        if charge_density is None:
            return None
        charge_density = {'charge_density': charge_density}
        return charge_density

    @property
    def magnetization_density(self):
        """
        Return the magnetization density.

        Returns
        -------
        magnetization_density : dict or ndarray
            If collinear spin calculations have been performed, a NumPy array containing
            the magnetization density in the unit cell in C order is returned. If however
            a non-collinear spin calculation have been performed, a dictionary is returned
            with keys `x`, `y` and `z`, each containing the same NumPy array, but for each
            direction.

        """
        magnetization_density = self._content_parser.magnetization_density
        if magnetization_density is None:
            return None
        magnetization_density = {'magnetization_density': magnetization_density}
        return magnetization_density
