"""
Common code for parsers.

------------------------
"""
from aiida.parsers.parser import Parser
from aiida.common.exceptions import NotExistent


class BaseParser(Parser):
    """Does common tasks all parsers carry out and provides convenience methods."""

    def __init__(self, node):
        super(BaseParser, self).__init__(node)
        self.out_folder = None

    def parse(self, **kwargs):
        """Set the out_folder attribute for use in extending parsers."""
        error_code = self.get_folder()
        if error_code is not None:
            return error_code
        return None

    def get_folder(self):
        """Convenient access to the retrieved folder."""
        try:
            _ = self.retrieved
            return None
        except NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

    def get_file(self, fname):
        """
        Convenient access to retrieved files.

        :param fname: name of the file
        :return: absolute path to the retrieved file
        """
        try:
            with self.retrieved.open(fname) as file_obj:
                ofname = file_obj.name
            return ofname
        except OSError:
            self.logger.warning(fname + ' not found in retrieved')
            return None
