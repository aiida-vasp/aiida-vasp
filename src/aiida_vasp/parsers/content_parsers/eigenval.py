"""
The ``EIGENVAL`` parser interface.

Contains the parsing interfaces to parsevasp used to parse ``EIGENVAL`` content.
"""

# pylint: disable=abstract-method
from typing import TextIO

from parsevasp.eigenval import Eigenval

from aiida_vasp.parsers.content_parsers.base import BaseFileParser


class EigenvalParser(BaseFileParser):
    """The parser interface that enables parsing of ``EIGENVAL`` content.

    The parser is triggered by using the ``eigenval-eigenvalues`` and/or ``eigenval-kpoints`` quantity key.
    The quantity keys ``eigenvalues`` and ``kpoints`` will on the other hand parse the
    eigenvalues and/or kpoints using the XML parser.

    """

    DEFAULT_SETTINGS = {'quantities_to_parse': ['eigenval-eigenvalues', 'eigenval-kpoints']}

    PARSABLE_QUANTITIES = {
        'eigenval-eigenvalues': {
            'inputs': [],
            'name': 'eigenvalues',
            'prerequisites': [],
        },
        'eigenval-kpoints': {
            'inputs': ['structure'],
            'name': 'kpoints',
            'prerequisites': ['structure'],
        },
    }

    def _init_from_handler(self, handler: TextIO) -> None:
        """Initialize a ``parsevasp`` object of ``Eigenval`` using a file like handler.

        :param handler: A file like object that provides the necessary ``EIGENVAL`` content to be parsed.
        :type handler: file-like object
        """

        try:
            self._content_parser = Eigenval(file_handler=handler, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    @property
    def eigenvalues(self):
        """Return the eigenvalue and metadata in the ``EIGENVAL``.

        :returns: A dict containing the keys ``eigenvalues`` and ``metadata``, which contain
                  a NumPy array of the eigenvalues and a dict with metadata, respectively.
        :rtype: dict
        """

        metadata = self._content_parser.get_metadata()
        _eigenvalues = self._content_parser.get_eigenvalues()
        eigenvalues = {'eigenvalues': _eigenvalues, 'metadata': metadata}

        return eigenvalues

    @property
    def kpoints(self):
        """Return the kpoints and metadata in the ``EIGENVAL``.

        :returns: A dict containing the keys ``kpoints`` and ``metadata``, which contain
                  a NumPy array of the k-points and a dict with metadata, respectively.
        :rtype: dict
        """

        metadata = self._content_parser.get_metadata()
        _kpoints = self._content_parser.get_kpoints()
        kpoints = {'kpoints': _kpoints, 'metadata': metadata}

        return kpoints
