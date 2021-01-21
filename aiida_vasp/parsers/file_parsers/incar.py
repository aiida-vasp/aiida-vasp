"""
INCAR parser.

-------------
The file parser that handles the parsing of INCAR files.
"""
from aiida.common import InputValidationError
from parsevasp.incar import Incar

from aiida_vasp.parsers.file_parsers.parser import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class


class IncarParser(BaseFileParser):
    """
    Parser for VASP INCAR format.

    This is a wrapper for the parsevasp.incar parser.

    The Parsing direction depends on whether the IncarParser is initialised with
    'path = ...' (read from file) or 'data = ...' (read from data).

    """

    PARSABLE_ITEMS = {
        'incar': {
            'inputs': [],
            'name': 'incar',
            'prerequisites': []
        },
    }

    def __init__(self, *args, **kwargs):
        super(IncarParser, self).__init__(*args, **kwargs)
        self.init_with_kwargs(**kwargs)

    def _init_with_data(self, data):
        """Initialize with a given AiiDA Dict instance."""
        if isinstance(data, get_data_class('dict')):
            self._data_obj = data
        else:
            self._logger.warning('Please supply an AiiDA Dict datatype for `data`.')
            self._data_obj = None
        self.parsable_items = self.__class__.PARSABLE_ITEMS
        self._parsed_data = {}

    @property
    def _parsed_object(self):
        """
        Return an instance of parsevasp.incar.Incar.

        Corresponds to the stored data in inputs.parameters.

        """

        incar_dict = self._data_obj.get_dict()

        try:
            return Incar(incar_dict=incar_dict, logger=self._logger)
        except SystemExit as error:
            raise InputValidationError(error.args[0])

    def _parse_file(self, inputs):
        """Create a DB Node from an INCAR file."""

        result = inputs
        result = {}

        if isinstance(self._data_obj, get_data_class('dict')):
            return {'incar': self._data_obj}

        try:
            incar = Incar(file_path=self._data_obj.path, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exitited abnormally. Returning None.')
            return {'incar': None}

        result['incar'] = incar.get_dict()
        return result

    @property
    def incar(self):
        return self.get_quantity('incar', {})['incar']
