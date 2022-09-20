"""
The ``KPOINTS`` parser interface.

-----------------------------
Contains the parsing interfaces to parsevasp used to parse ``KPOINTS`` content.
"""
import numpy as np

from parsevasp.kpoints import Kpoints, Kpoint
from aiida_vasp.parsers.content_parsers.base import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class


class KpointsParser(BaseFileParser):
    """The parser interface that enables parsing of ``KPOINTS`` content.

    The parser is triggered by using the ``kpoints-kpoints`` quantity key. The quantity key ``kpoints``
    will on the other hand parse the k-points using the XML parser.

    """

    DEFAULT_SETTINGS = {'quantities_to_parse': ['kpoints-kpoints']}

    PARSABLE_QUANTITIES = {
        'kpoints-kpoints': {
            'inputs': [],
            'name': 'kpoints',
            'prerequisites': [],
        },
    }

    def _init_from_handler(self, handler):
        """Initialize using a file like handler.

        Parameters
        ----------
        handler : object
            A file like object that provides the necessary content to be parsed.

        """

        try:
            self._content_parser = Kpoints(file_handler=handler, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    def _init_from_data(self, data):
        """Initialize using AiiDA KpointsData."""

        if isinstance(data, get_data_class('array.kpoints')):
            self._content_data = data
        else:
            raise TypeError('The supplied AiiDA data structure is not a KpointsData.')

    @property
    def kpoints(self):
        """
        Return kpoints that is ready to be consumed by the the AiiDA ``KpointsData``.

        AiiDA does not support the line mode used in VASP, so we give a warning that parsing
        this is not supported.

        Returns
        -------
        aiida_kpoints : dict
            A dict that contain keys ``comment``, ``divisions``, ``shifts``, ``points``, ``tetra``,
            ``tetra_volume``, ``mode`` ``centering``, ``num_kpoints``, ``weights`` and ``cartesian``
            which are compatible with consumption of the initialization of the AiiDA KpointsData.

        """

        aiida_kpoints = parsevasp_to_aiida(self._content_parser, self._logger)

        return aiida_kpoints

    def _content_data_to_content_parser(self):
        """
        Convert an AiiDA ``KpointsData`` to a content parser instance of ``Kpoints`` from ``parsevasp``.

        Returns
        -------
        content_parser : object
            An instance of ``Kpoints`` from ``parsevasp``.

        """
        try:
            # Check if the ``KpointsData`` contain a mesh.
            _ = self._content_data.get_attribute('mesh')
            mode = 'automatic'
        except AttributeError:
            pass

        try:
            # Check to see if the ``KpointsData`` contain an explicit k-point list.
            _ = self._content_data.get_attribute('array|kpoints')
            mode = 'explicit'
        except AttributeError:
            pass

        kpoints_dict = {}
        for keyword in [
                'comment', 'divisions', 'shifts', 'points', 'tetra', 'tetra_volume', 'mode', 'centering', 'num_kpoints',
                'generating_vectors'
        ]:
            kpoints_dict[keyword] = None

        kpoints_dict.update(getattr(self, '_get_kpointsdict_' + mode)(self._content_data))

        # We brake hard if ``parsevasp`` fail here. If we can not write we will not try another parser.
        content_parser = Kpoints(kpoints_dict=kpoints_dict, logger=self._logger)

        return content_parser

    def _get_kpointsdict_explicit(self, kpoints_data):  # pylint: disable=no-self-use
        """
        Turn Aiida ``KpointData`` into a k-points dictionary with explicit generation of points.

        Parameters
        ----------
        kpoints_data : object
            An AiiDA ``KpointsData`` object containing explicit k-point sets.

        Returns
        -------
        kpoints_dict : dict
            A dictionary that can be used to initialize a ``parsevasp`` ``Kpoints`` instance.

        """
        kpoints_dict = {}

        kpts = []
        try:
            points, weights = kpoints_data.get_kpoints(also_weights=True)
        except AttributeError:
            points = kpoints_data.get_kpoints()
            weights = None
        for index, point in enumerate(points):
            if weights is not None:
                kpt = Kpoint(point, weight=weights[index])
            else:
                # No weights supplied, so set them to 1.0
                kpt = Kpoint(point, weight=1.0)
            kpts.append(kpt)
        kpoints_dict['points'] = kpts
        kpoints_dict['mode'] = 'explicit'
        kpoints_dict['num_kpoints'] = len(kpts)

        return kpoints_dict

    @staticmethod
    def _get_kpointsdict_automatic(kpointsdata):
        """
        Turn Aiida ``KpointData`` into a k-point dictionary with automatic generation of points.

        Parameters
        ----------
        kpoints_data : object
            An AiiDA ``KpointsData`` object containing meshed k-point sets.

        Returns
        -------
        kpoints_dict : dict
            A dictionary that can be used to initialize a ``parsevasp`` ``Kpoints`` instance.

        """

        kpoints_dict = {}
        # Automatic mode
        mesh = kpointsdata.get_kpoints_mesh()
        kpoints_dict['divisions'] = mesh[0]
        kpoints_dict['shifts'] = mesh[1]
        kpoints_dict['mode'] = 'automatic'
        # Here we need to make a choice, so should add more to AiiDA to make this better defined
        kpoints_dict['centering'] = 'Gamma'
        kpoints_dict['num_kpoints'] = 0

        return kpoints_dict


def parsevasp_to_aiida(kpoints, logger):
    """``parsevasp`` to AiiDA conversion.

    Generate an AiiDA data structure that can be consumed by ``KpointsData`` on initialization
    from the ``parsevasp`` instance of the ``Kpoints`` class.

    Parameters
    ----------
    kpoints : object
        An instance of the ``Kpoints`` class in ``parsevasp``.

    Returns
    -------
    kpoints_dict : dict
        A dictionary representation which are ready to be used when creating an
        AiiDA ``KpointsData`` instance.

    """

    if kpoints.entries.get('mode') == 'line':
        # AiiDA does not support line mode
        logger.warning('The read KPOINTS contained line mode which is' 'not supported. Returning None.')
        return None

    # Fetch a dictionary containing the k-points information
    kpoints_dict = kpoints.get_dict()

    # Now unpack points, weights and check direct versus cartesian. Set to
    # None if mode is automatic.
    points = []
    weights = []
    cartesian = []
    for key, value in kpoints_dict.items():
        if key == 'points':
            if value is not None:
                for item in value:
                    points.append(item[0])
                    weights.append(item[1])
                    # AiiDA wants cartesian and not direct flags, so revert
                    cartesian.append(not item[2])
            else:
                points = None
                weights = None
                cartesian = None

    # Make sure weights is ndarray
    if weights is not None:
        weights = np.array(weights)

    # Check that we only have similar elements in the direct list as
    # AiiDA can only work with all points being either in direct or cartesian
    # coordinates.
    if cartesian is not None:
        if not cartesian.count(cartesian[0]) == len(cartesian):
            raise ValueError('Different coordinate systems have been detected among the k-points.')
        cartesian = cartesian[0]

    # Modify dict to AiiDA spec
    kpoints_dict['points'] = points
    kpoints_dict['weights'] = weights
    kpoints_dict['cartesian'] = cartesian

    return kpoints_dict
