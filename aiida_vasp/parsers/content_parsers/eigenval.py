"""
The ``EIGENVAL`` parser interface.

----------------------------------
Contains the parsing interfaces to parsevasp used to parse ``EIGENVAL`` content.
"""
# pylint: disable=abstract-method
from aiida_vasp.parsers.content_parsers.base import BaseFileParser

from parsevasp.eigenval import Eigenval


class EigenvalParser(BaseFileParser):
    """The parser interface that enables parsing of ``EIGENVAL`` content.

    The parser is triggered by using the ``eigenval-bands`` and/or ``eigenval-kpoints`` quantity key.
    The quantity keys ``bands`` and ``kpoints`` will on the other hand parse the
    eigenvalues and/or kpoints using the XML parser.

    """

    DEFAULT_OPTIONS = {'quantities_to_parse': ['eigenval-bands', 'eigenval-kpoints']}

    PARSABLE_QUANTITIES = {
        'eigenval-bands': {
            'inputs': [],
            'name': 'bands',
            'prerequisites': [],
        },
        'eigenval-kpoints': {
            'inputs': ['structure'],
            'name': 'kpoints',
            'prerequisites': ['structure'],
        },
    }

    def _init_from_handler(self, handler):
        """Initialize a ``parsevasp`` object of ``Eigenval`` using a file like handler.

        Parameters
        ----------
        handler : object
            A file like object that provides the necessary ``EIGENVAL`` content to be parsed.

        """

        try:
            self._content_parser = Eigenval(file_handler=handler, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    @property
    def bands(self):
        """Return the eigenvalue and metadata in the ``EIGENVAL``.

        Returns
        -------
        bands : dict
            A dict containing the keys ``bands`` and ``metadata``, which contain
            a NumPy array of the eigenvalues and a dict with metadata, respectively.

        """

        metadata = self._content_parser.get_metadata()
        _bands = self._content_parser.get_bands()
        bands = {'bands': _bands, 'metadata': metadata}

        return bands

    @property
    def kpoints(self):
        """Return the kpoints and metadata in the ``EIGENVAL``.

        Returns
        -------
        kpoints : dict
            A dict containing the keys ``kpoints`` and ``metadata``, which contain
            a NumPy array of the k-points and a dict with metadata, respectively.

        """

        metadata = self._content_parser.get_metadata()
        _kpoints = self._content_parser.get_kpoints()
        kpoints = {'kpoints': _kpoints, 'metadata': metadata}

        return kpoints
