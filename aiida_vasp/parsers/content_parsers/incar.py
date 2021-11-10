"""
INCAR parser.

-------------
Contains the parsing interfaces to parsevasp used to parse INCAR.
"""

from aiida_vasp.parsers.content_parsers.parser import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class

from parsevasp.incar import Incar


class IncarParser(BaseFileParser):
    """The parser interface that enables parsing of POSCAR/CONTCAR.

    The parser is triggered by using the `incar` quantity key.

    """

    DEFAULT_OPTIONS = {'quantities_to_parse': ['incar']}

    PARSABLE_QUANTITIES = {
        'incar': {
            'inputs': [],
            'name': 'incar',
            'prerequisites': []
        },
    }

    def _init_from_handler(self, handler):
        """Initialize using a file like handler."""

        try:
            self._content_parser = Incar(file_handler=handler, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    def _init_from_data(self, data):
        """Initialize using an AiiDA Dict data type."""

        if isinstance(data, get_data_class('dict')):
            self._content_data = data
        else:
            raise TypeError('The supplied AiiDA data structure is not a Dict.')

    @property
    def incar(self):
        """Return the parameters in the INCAR.

        Returns
        -------
        params : dict or None
            A dictionary containing the parameter tags as keys and its settings as values.
            None is returned if the quantity can not be parsed.

        """
        if self._content_parser is not None:
            params = self._content_parser.get_dict()
            return params
        return None

    def _content_data_to_content_parser(self):
        """
        Convert an AiiDA Dict to a content parser instance of Incar from parsevasp.

        Returns
        -------
        incar : object
            An instance of Incar from parsevasp.

        """
        dictionary = self._content_data.get_dict()

        # We brake hard if parsevasp fail here. If we can not write we will not try another parser.
        incar = Incar(incar_dict=dictionary, logger=self._logger)

        return incar
