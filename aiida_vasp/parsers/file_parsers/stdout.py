"""
STDOUT parser.

--------------
File parser that handles parsing of errors found in STDOUT
"""

from aiida_vasp.parsers.file_parsers.parser import BaseFileParser, SingleFile
from aiida_vasp.parsers.error_handling import ErrorScanner

DEAULT_OPTIONS = {'quantities_to_parse': ['stdout_error']}


class StdoutParser(BaseFileParser):
    """Parser used for parsing errors from STDOUT"""
    PARSABLE_ITEMS = {
        'stdout_error': {
            'inputs': [],
            'name': 'stdout_error',
            'prerequisites': [],
        }
    }

    def __init__(self, *args, **kwargs):
        super(StdoutParser, self).__init__(*args, **kwargs)
        self.scanner = None
        self.init_with_kwargs(**kwargs)

    def _init_with_file_path(self, path):
        """Init with a file path"""
        self._parsed_data = {}
        self.parsable_items = self.__class__.PARSABLE_ITEMS
        self._data_obj = SingleFile(path=path)

    def _init_with_data(self, data):
        """Init with SingleFileData"""
        self.parsable_items = self.__class__.PARSABLE_ITEMS
        self._init_with_file_path(data.get_file_abs_path())

    def _parse_file(self, inputs):
        """Parse the STDOUT file"""

        del inputs
        quantities_to_parse = DEAULT_OPTIONS.get('quantities_to_parse')
        if self.settings is not None:
            quantities_to_parse = self.settings.quantities_to_parse

        result = {}
        if 'stdout_error' in quantities_to_parse:
            with open(self._data_obj.path, 'r') as fhandle:
                scanner = ErrorScanner(fhandle, ftype='STDOUT')
                scanner.scan()
                errors = scanner.get_errors()
            # Construct the error dict from ErrorRecord named tuple
            error_dict_keys = ['shortname', 'suggestion', 'message', 'critical']
            error_dicts = []
            for error in errors:
                error_dict = {key: getattr(error, key) for key in error_dict_keys}
                if error_dict['suggestion'] is None:
                    # Remove the suggestion if it has no content
                    error_dict.pop('suggestion')
                error_dicts.append(error_dict)
            result['stdout_error'] = error_dicts

        return result
