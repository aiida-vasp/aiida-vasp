"""
KPOINTS parser.

---------------
Contains the parsing interfaces to parsevasp used to parse KPOINTS.
"""

from parsevasp.kpoints import Kpoints, Kpoint
from aiida_vasp.parsers.content_parsers.parser import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class


class KpointsParser(BaseFileParser):
    """The parser interface that enables parsing of KPOINTS.

    The parser is triggered by using the `kpoints-kpoints` quantity key. The quantity key `kpoints`
    will on the other hand parse the k-points using the XML parser.

    """

    DEFAULT_OPTIONS = {'quantities_to_parse': ['kpoints-kpoints']}

    PARSABLE_QUANTITIES = {
        'kpoints-kpoints': {
            'inputs': [],
            'name': 'kpoints',
            'prerequisites': [],
        },
    }

    def _init_from_handler(self, handler):
        """Initialize using a file like handler."""

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
        Return kpoints that is ready to be consumed by the the AiiDA KpointsData.

        AiiDA does not support the line mode used in VASP, so we give a warning that parsing
        this is not supported.

        Returns
        -------
        aiida_kpoints : dict
            A dict that contain keys `comment`, `divisions`, `shifts`, `points`, `tetra`,
            `tetra_volume`, `mode` `centering` and `num_kpoints`, which are compatible
            with consumption of the initialization of the AiiDA StructureData.

        """

        aiida_kpoints = None
        if self._content_parser.entries.get('mode') == 'line':
            self._logger.warning('The read KPOINTS contained line mode which is' 'not supported. Returning None.')
        else:
            aiida_kpoints = self._content_parser.get_dict()

        return aiida_kpoints

    def _content_data_to_content_parser(self):
        """
        Convert an AiiDA KpointsData to a content parser instance of Kpoints from parsevasp.

        Returns
        -------
        kpoints : object
            An instance of Kpoints from parsevasp.

        """
        try:
            # Check if the KpointsData contain a mesh.
            _ = self._content_data.get_attribute('mesh')
            mode = 'automatic'
        except AttributeError:
            pass

        try:
            # Check to see if the KpointsData contain an explicit k-point list.
            _ = self._content_data.get_attribute('array|kpoints')
            mode = 'explicit'
        except AttributeError:
            pass

        kpoints_dict = {}
        for keyword in ['comment', 'divisions', 'shifts', 'points', 'tetra', 'tetra_volume', 'mode', 'centering', 'num_kpoints']:
            kpoints_dict[keyword] = None

            kpoints_dict.update(getattr(self, '_get_kpointsdict_' + mode)(self._content_data))

        # We brake hard if parsevasp fail here. If we can not write we will not try another parser.
        kpoints = Kpoints(kpoints_dict=kpoints_dict, logger=self._logger)

        return kpoints

    def _get_kpointsdict_explicit(self, kpoints_data):
        """
        Turn Aiida KpointData into an 'explicit' kpoints dictionary.

        Parameters
        ----------
        kpoints_data : object
            An AiiDA KpointsData object containing explicit k-point sets.

        Returns
        -------
        kpoints_dict : dict
            A dictionary that can be used to initialize a parsevasp Kpoints instance.

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
                kpt = Kpoint(point, weight=weights[index], logger=self._logger)
            else:
                # no weights supplied, so set them to 1.0
                kpt = Kpoint(point, weight=1.0, logger=self._logger)
            kpts.append(kpt)
        kpoints_dict['points'] = kpts
        kpoints_dict['mode'] = 'explicit'
        kpoints_dict['num_kpoints'] = len(kpts)

        return kpoints_dict

    @staticmethod
    def _get_kpointsdict_automatic(kpointsdata):
        """
        Turn Aiida KpointData into an 'automatic' kpoints dictionary.

        Parameters
        ----------
        kpoints_data : object
            An AiiDA KpointsData object containing meshed k-point sets.

        Returns
        -------
        kpoints_dict : dict
            A dictionary that can be used to initialize a parsevasp Kpoints instance.

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

    # @staticmethod
    # def _get_kpointsdata_explicit(kpoints_dict):
    #     """Turn an 'explicit' kpoints dictionary into Aiida KpointsData."""
    #     kpout = get_data_class('array.kpoints')()

    #     kpoints = kpoints_dict.get('points')
    #     cartesian = not kpoints[0].get_direct()
    #     kpoint_list = []
    #     weights = []
    #     for kpoint in kpoints:
    #         kpoint_list.append(kpoint.get_point().tolist())
    #         weights.append(kpoint.get_weight())

    #     if weights[0] is None:
    #         weights = None

    #     kpout.set_kpoints(kpoint_list, weights=weights, cartesian=cartesian)

    #     return kpout

    # @staticmethod
    # def _get_kpointsdata_automatic(kpoints_dict):
    #     """Turn an 'automatic' kpoints dictionary into Aiida KpointsData."""
    #     kpout = get_data_class('array.kpoints')()

    #     mesh = kpoints_dict.get('divisions')
    #     shifts = kpoints_dict.get('shifts')
    #     kpout.set_kpoints_mesh(mesh, offset=shifts)

    #     return kpout
