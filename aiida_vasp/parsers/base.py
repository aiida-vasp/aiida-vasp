"""
Common code for parsers.

------------------------
"""
import os
from contextlib import contextmanager

from aiida.parsers.parser import Parser
from aiida.common.exceptions import NotExistent
from aiida.repository.common import FileType


class BaseParser(Parser):
    """Does common tasks all parsers carry out and provides convenience methods."""

    def __init__(self, node):
        super().__init__(node)
        self._retrieved_content = None
        self._retrieved_temporary = None

    def parse(self, **kwargs):
        """Check the folders and set the retrieved_content for use in extending parsers."""
        error_code = self._compose_retrieved_content(kwargs)
        if error_code is not None:
            return error_code
        return None

    def _compose_retrieved_content(self, parser_kwargs=None):  # pylint:disable=too-many-branches
        """
        Convenient check to see if the retrieved and retrieved temporary folder is present.

        This routine also builds a dictionary containing the content of both the retrieved folder and
        the retrieved_temporary folder, accessible from retrieved_content. The error for the temporary
        folder takes presence as this is the one we mostly rely on.
        """
        retrieved = {}
        exit_code_permanent = self._check_folder()
        if exit_code_permanent is None:
            # Retrieved folder exists, add content and tag to dictionary
            for retrieved_object in self.retrieved.list_objects():
                # Add sub directory objects - this only treat for one extra level
                if retrieved_object.file_type == FileType.DIRECTORY:
                    for subobj in self.retrieved.list_objects(retrieved_object.name):
                        retrieved[retrieved_object.name + '/' + subobj.name] = {'path': '', 'status': 'permanent'}
                else:
                    retrieved[retrieved_object.name] = {'path': '', 'status': 'permanent'}

        exit_code_temporary = None
        if parser_kwargs is not None:
            exit_code_temporary = self._set_retrieved_temporary(parser_kwargs)
            if exit_code_temporary is None:
                # Retrieved_temporary folder exists, add content and tag to dictionary
                # Notice https://github.com/aiidateam/aiida-core/issues/3502. That is why
                # we need to store the path and use listdir etc.
                for retrieved_object in os.listdir(self._retrieved_temporary):
                    abspath = os.path.join(self._retrieved_temporary, retrieved_object)
                    # Add sub directory objects - this only treat for one extra level
                    if os.path.isdir(abspath):
                        for subobj in os.listdir(abspath):
                            retrieved[os.path.join(retrieved_object, subobj)] = {'path': self._retrieved_temporary, 'status': 'temporary'}
                    else:
                        retrieved[retrieved_object] = {'path': self._retrieved_temporary, 'status': 'temporary'}

        # Check if there are other objects than the AiiDA generated scheduler objects in retrieved and
        # if there are any objects in the retrieved_temporary. If not, return an error.
        aiida_required_objects = [self.node.get_attribute('scheduler_stderr'), self.node.get_attribute('scheduler_stdout')]
        # Check if have some missing objects that we require to be present.
        vasp_output_objects_present = False
        for name in retrieved:
            if name not in aiida_required_objects:
                vasp_output_objects_present = True
        if not vasp_output_objects_present:
            return self.exit_codes.ERROR_VASP_DID_NOT_EXECUTE

        # Store the retrieved content.
        self._retrieved_content = retrieved
        # At least one of the folders should be present.
        if exit_code_permanent is None or exit_code_temporary is None:
            return None
        # Both are not present, exit code of the temporary folder take precedence as this is
        # the one we typically use the most.
        return exit_code_temporary if not None else exit_code_permanent

    def _set_retrieved_temporary(self, parser_kwargs):
        """
        Check the presence of retrieved_temporary folder.

        This folder is typically set as a SandboxFolder if there is anything in the retrieve_temporary_list
        and deleted when the parsing is completed.
        """
        self._retrieved_temporary = parser_kwargs.get('retrieved_temporary_folder', None)
        if self._retrieved_temporary is None:
            return self.exit_codes.ERROR_NO_RETRIEVED_TEMPORARY_FOLDER
        return None

    def _check_folder(self):
        """Check the presence of the retrieved folder."""
        try:
            _ = self.retrieved
            return None
        except NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

    @contextmanager
    def _get_handler(self, name, mode, encoding=None):
        """
        Access the handler of retrieved and retrieved_temporary objects. This is passed
        down to the parers where the content is analyzed.

        Since their interface is presently different, we created this wrapper to
        make it uniform. See also https://github.com/aiidateam/aiida-core/issues/3502.

        Furthermore, since we now use context managers for e.g. file handler we needed
        to also make this method a contextmanager. By using the decorator we do not have
        to make it a class and define the __enter__ and __exit__ methods.

        We also allow opening the file in different modes.

        :param name: name of the object
        :param mode: the open mode to use for the respective parser, typically 'r' or 'rb'.
        :param encoding: the encoding to be used, if binary mode, this is ignored, otherwise defaults to utf8 if not given.
        :returns: a yielded handler for the object
        :rtype: object
        """

        if 'b' not in mode and encoding is None:
            # If we do not use binary, set default is encoding not proided
            encoding = 'utf8'
        try:
            if self._retrieved_content[name]['status'] == 'permanent':
                # For the permanent content to be parsed we can use the fact that
                # self.retrieved is a FolderData datatype in AiiDA.
                try:
                    with self.retrieved.open(name, mode) as handler:  # pylint: disable=unspecified-encoding
                        yield handler
                except OSError:
                    self.logger.warning(name + ' not found in retrieved')
                    yield None
            else:
                # For the temporary content to be parsed we have to use the regular
                # folder approach for now.
                # See https://github.com/aiidateam/aiida-core/issues/3502.
                path = os.path.join(self._retrieved_content[name]['path'], name)
                try:
                    with open(path, mode, encoding=encoding) as handler:
                        yield handler
                except OSError:
                    self.logger.warning(name + ' not found in retrieved_temporary')
                    yield None
        except KeyError:
            yield None


def list_files_recursive(retrieved, top_level=''):
    """
    Recursively list the content of a retrieved FolderData node.
    """
    object_paths = []
    for obj in retrieved.list_objects(top_level):
        if obj.file_type == FileType.FILE:
            object_paths.append(os.path.join(top_level, obj.name))
        elif obj.file_type == FileType.DIRECTORY:
            object_paths.extend(
                [os.path.join(top_level, path) for path in list_files_recursive(retrieved, os.path.join(top_level, obj.name))])

    return object_paths
