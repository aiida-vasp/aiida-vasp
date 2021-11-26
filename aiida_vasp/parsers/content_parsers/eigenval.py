"""
EIGENVAL parser.

----------------
Contains the parsing interfaces to parsevasp used to parse EIGENVAL.
"""
from aiida_vasp.parsers.content_parsers.base import BaseFileParser

from parsevasp.eigenval import Eigenval


class EigenvalParser(BaseFileParser):
    """The parser interface that enables parsing of EIGENVAL.

    The parser is triggered by using the `eigenval-bands` and/or `eigenval-kpoints` quantity key.
    The quantity keys `bands` and `kpoints` will on the other hand parse the
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
        """Initialize using a file like handler."""

        try:
            self._content_parser = Eigenval(file_handler=handler, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    @property
    def bands(self):
        """
        Return the eigenvalue and metadata in the EIGENVAL.

        Returns
        -------
        bands : dict
            A dict containing the tkeys `bands` and `metadata`, which contain
            a NumPy array and a dict, respectively.

        """
        metadata = self._content_parser.get_metadata()
        _bands = self._content_parser.get_bands()
        bands = {'bands': _bands, 'metadata': metadata}

        return bands

    @property
    def kpoints(self):
        metadata = self._content_parser.get_metadata()
        _kpoints = self._content_parser.get_kpoints()
        kpoints = {'kpoints': _kpoints, 'metadata': metadata}

        return kpoints
