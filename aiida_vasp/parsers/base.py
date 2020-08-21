"""
Common code for parsers.

------------------------
"""
import os

from aiida.parsers.parser import Parser
from aiida.common.exceptions import NotExistent


class BaseParser(Parser):
    """Does common tasks all parsers carry out and provides convenience methods."""

    def __init__(self, node):
        super(BaseParser, self).__init__(node)
        self.retrieved_content = None
        self.retrieved_temporary = None

    def parse(self, **kwargs):
        """Check the folders and set the retrieved_content for use in extending parsers."""
        error_code = self.check_folders(kwargs)
        if error_code is not None:
            return error_code
        return None

    def check_folders(self, parser_kwargs=None):
        """
        Convenient check to see if the retrieved and retrieved temp folder is present.

        This routine also builds a dictionary containing the content of both the retrieved folder and
        the retrieved_temporary folder, accessible from retrieved_content. The error for the temporary
        folder takes presence as this is the one we mostly rely on.
        """

        retrieved = {}
        exit_code_permanent = self.check_folder()
        if exit_code_permanent is None:
            # Retrieved folder exists, add content and tag to dictionary
            for retrieved_file in self.retrieved.list_objects():
                retrieved[retrieved_file.name] = {'path': '', 'status': 'permanent'}

        exit_code_temporary = None
        if parser_kwargs is not None:
            exit_code_temporary = self.check_temporary_folder(parser_kwargs)
            if exit_code_temporary is None:
                # Retrieved_temporary folder exists, add content and tag to dictionary
                for retrieved_file in os.listdir(self.retrieved_temporary):
                    retrieved[retrieved_file] = {'path': self.retrieved_temporary, 'status': 'temporary'}

        # Check if there are other files than the AiiDA generated scheduler files in retrieved and
        # if there are any files in the retrieved_temporary. If not, return an error.
        aiida_required_files = [self.node.get_attribute('scheduler_stderr'), self.node.get_attribute('scheduler_stdout')]
        vasp_output_files_present = False
        for file_name in retrieved:
            if file_name not in aiida_required_files:
                vasp_output_files_present = True

        if not vasp_output_files_present:
            return self.exit_codes.ERROR_VASP_DID_NOT_EXECUTE

        # Store the retrieved content
        self.retrieved_content = retrieved
        # OK if a least one of the folders are present
        if exit_code_permanent is None or exit_code_temporary is None:
            return None
        # Both are not present, exit code of the temporary folder take precedence
        return exit_code_temporary if not None else exit_code_permanent

    def check_temporary_folder(self, parser_kwargs):
        """Convenient check of the retrieved_temporary folder."""
        self.retrieved_temporary = parser_kwargs.get('retrieved_temporary_folder', None)
        if self.retrieved_temporary is None:
            return self.exit_codes.ERROR_NO_RETRIEVED_TEMPORARY_FOLDER
        return None

    def check_folder(self):
        """Convenient check of the retrieved folder."""
        try:
            _ = self.retrieved
            return None
        except NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

    def get_file(self, fname):
        """
        Convenient access to retrieved and retrieved_temporary files.

        :param fname: name of the file
        :return: absolute path to the retrieved file
        """

        try:
            if self.retrieved_content[fname]['status'] == 'permanent':
                try:
                    with self.retrieved.open(fname) as file_obj:
                        ofname = file_obj.name
                    return ofname
                except OSError:
                    self.logger.warning(fname + ' not found in retrieved')
                    return None
            else:
                path = self.retrieved_content[fname]['path']
                file_path = os.path.join(path, fname)
                try:
                    with open(file_path, 'r') as file_obj:
                        ofname = file_path
                    return ofname
                except OSError:
                    self.logger.warning(fname + ' not found in retrieved_temporary')
                    return None
        except KeyError:
            return None
